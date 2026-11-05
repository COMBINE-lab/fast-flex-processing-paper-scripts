"""
Shared I/O for the Flex concordance analysis: load one (dataset, sample) count matrix
from each tool onto a common cell/gene axis.

Axis conventions
----------------
cells : the 16 nt cell barcode. Within one Flex *sample* the probe barcode is fixed, so
        16 nt is a unique key. Cell Ranger writes cell-BC+probe-BC (24 nt for Flex v1,
        26 nt for v2); simpleaf writes ``<sample>_<16nt>``; cyto writes ``<16nt>-<sample>``.
genes : Ensembl gene_id. cyto reports gene *names*, harmonized through the probe-set CSV.

Cell Ranger emits one file per sample except when the multi config declares a single
pooled sample spanning several probe barcodes (320k_GEMX: 16 probe barcodes in one
sample, 289,331 cells over only 127,235 distinct 16 nt cores). There the per-sample cell
set is recovered by splitting on the probe-barcode suffix, which is column 2 of the
sample-barcode list (Cell Ranger reports the translated form).
"""
import glob
import gzip
import json
import os
import re
import shutil
import sys
import time

import numpy as np
import scipy.io
import scipy.sparse as sp

BC16 = re.compile(r"[ACGTN]{16}")
GEX = "Gene Expression"


def _open(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def norm_bc(s):
    """-> the 16 nt cell-barcode core of any tool's barcode label."""
    m = BC16.search(s.strip())
    return m.group(0) if m else s.strip()


def _find(d, *names):
    for n in names:
        hits = sorted(glob.glob(os.path.join(d, n)))
        if hits:
            return hits[0]
    return None


# ---------------------------------------------------------------- loaders

def load_10x_dir(mtx_dir, gene_col=0):
    """10x mtx triplet -> (csr cells×genes, barcodes[16nt], gene labels)."""
    mtx = _find(mtx_dir, "matrix.mtx.gz", "matrix.mtx")
    bcs = _find(mtx_dir, "barcodes.tsv.gz", "barcodes.tsv")
    fts = _find(mtx_dir, "features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv")
    if not (mtx and bcs and fts):
        return None
    with _open(mtx) as fh:
        M = scipy.io.mmread(fh).tocsr()           # 10x layout: genes × cells
    barcodes = [norm_bc(l) for l in _open(bcs).read().splitlines()]
    lines = [l.split("\t") for l in _open(fts).read().splitlines()]
    genes = [p[gene_col] for p in lines]
    ftypes = [p[2] if len(p) > 2 else GEX for p in lines]
    M = M.T.tocsr()                               # -> cells × genes
    keep = [i for i, t in enumerate(ftypes) if t == GEX]
    if len(keep) != len(genes):
        M = M[:, keep].tocsr()
        genes = [genes[i] for i in keep]
    return M, barcodes, genes


def load_10x_h5(path):
    """Cell Ranger .h5 -> (csr cells×genes, raw barcode strings, gene_ids). GEX only."""
    import h5py

    with h5py.File(path, "r") as f:
        g = f["matrix"]
        data, indices, indptr = g["data"][:], g["indices"][:], g["indptr"][:]
        n_genes, n_cells = g["shape"][:]
        M = sp.csc_matrix((data, indices, indptr), shape=(n_genes, n_cells)).T.tocsr()
        barcodes = [b.decode() for b in g["barcodes"][:]]
        ids = [x.decode() for x in g["features"]["id"][:]]
        ftypes = [x.decode() for x in g["features"]["feature_type"][:]]
    keep = [i for i, t in enumerate(ftypes) if t == GEX]
    if len(keep) != len(ids):
        M = M[:, keep].tocsr()
        ids = [ids[i] for i in keep]
    return M, barcodes, ids


def load_cellranger(path, probe_suffix=None):
    """Load a Cell Ranger filtered matrix (.h5 file or mtx dir).

    ``probe_suffix`` keeps only the cells whose barcode carries that probe barcode —
    needed when one Cell Ranger sample pools several probe barcodes.
    """
    if os.path.isdir(path):
        sub = os.path.join(path, "sample_filtered_feature_bc_matrix")
        loaded = load_10x_dir(sub if os.path.isdir(sub) else path)
        if loaded is None:
            h5 = _find(path, "*sample_filtered_feature_bc_matrix.h5")
            loaded = load_10x_h5(h5) if h5 else None
        if loaded is None:
            return None
        M, raw_bc, genes = loaded
    else:
        M, raw_bc, genes = load_10x_h5(path)

    cores = [b.split("-")[0] for b in raw_bc]
    if probe_suffix:
        rows = [i for i, c in enumerate(cores) if c[16:] == probe_suffix]
        if not rows:
            return None
        M = M[rows, :].tocsr()
        cores = [cores[i] for i in rows]
    return M, [c[:16] for c in cores], genes


SPLIT_ABOVE_BYTES = 2 << 30      # stream+cache instead of loading whole .mtx


def _split_alevin_by_sample(mtx, row_lbl, n_genes, cache_dir, chunk_rows=20_000_000):
    """One bounded-memory pass over a large quants_mat.mtx -> one .npz per sample.

    A multiplexed quant of a 16-plex library carries millions of unfiltered barcodes
    (320k_GEMX: 4.9 M rows, tens of GB of MatrixMarket text, ~10^9 entries). The entries
    are read in chunks and appended to per-sample scratch files, so peak memory is one
    chunk (~320 MB) rather than the whole matrix, and the parse happens once instead of
    once per sample and analysis stage.
    """
    import pandas as pd

    os.makedirs(cache_dir, exist_ok=True)
    order, index_of = [], {}
    sample_of = np.empty(len(row_lbl), dtype=np.int32)
    local_row = np.empty(len(row_lbl), dtype=np.int64)
    counts = []
    for gi, lbl in enumerate(row_lbl):
        smp = lbl.rsplit("_", 1)[0]
        si = index_of.get(smp)
        if si is None:
            si = index_of[smp] = len(order)
            order.append(smp)
            counts.append(0)
        sample_of[gi] = si
        local_row[gi] = counts[si]
        counts[si] += 1

    scratch = [(os.path.join(cache_dir, f"_{si}.rows"),
                os.path.join(cache_dir, f"_{si}.cols"),
                os.path.join(cache_dir, f"_{si}.vals")) for si in range(len(order))]
    for paths in scratch:
        for f in paths:
            open(f, "wb").close()

    with _open(mtx) as fh:
        for line in fh:                       # skip banner/comments, then the dims line
            if not line.startswith("%"):
                break
        reader = pd.read_csv(fh, sep=r"\s+", header=None, chunksize=chunk_rows,
                             names=("r", "c", "v"),
                             dtype={"r": np.int64, "c": np.int32, "v": np.float32})
        for chunk in reader:
            gi = chunk["r"].to_numpy() - 1
            si = sample_of[gi]
            lr = local_row[gi]
            cols = chunk["c"].to_numpy() - 1
            vals = chunk["v"].to_numpy()
            for s in np.unique(si):
                m = si == s
                rp, cp, vp = scratch[s]
                with open(rp, "ab") as f:
                    lr[m].astype(np.int64).tofile(f)
                with open(cp, "ab") as f:
                    cols[m].astype(np.int32).tofile(f)
                with open(vp, "ab") as f:
                    vals[m].astype(np.float32).tofile(f)

    for si, smp in enumerate(order):
        rp, cp, vp = scratch[si]
        M = sp.coo_matrix(
            (np.fromfile(vp, dtype=np.float32),
             (np.fromfile(rp, dtype=np.int64), np.fromfile(cp, dtype=np.int32))),
            shape=(counts[si], n_genes)).tocsr()
        sp.save_npz(os.path.join(cache_dir, smp + ".npz"), M)
        for f in (rp, cp, vp):
            os.remove(f)
    with open(os.path.join(cache_dir, "_done"), "w") as f:
        f.write("\n".join(order) + "\n")


def _ensure_split(mtx, alv, row_lbl, n_genes, timeout=7200):
    """Build the per-sample cache exactly once, even with many readers in parallel.

    The driver runs several samples of one dataset concurrently, so without this they would
    all parse the same multi-GB matrix at the same time and interleave writes into one cache
    directory. One process wins the lock and builds into a private temp dir that is renamed
    into place atomically; the others wait for the marker.
    """
    cache = os.path.join(alv, "_by_sample")
    done = os.path.join(cache, "_done")
    if os.path.exists(done):
        return cache
    lock = os.path.join(alv, "_by_sample.lock")
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        waited = 0
        while not os.path.exists(done):        # another process is building it
            time.sleep(5)
            waited += 5
            if waited > timeout:
                raise TimeoutError(f"waited {waited}s for {done}")
        return cache
    try:
        os.write(fd, f"{os.getpid()}\n".encode())
        os.close(fd)
        tmp = cache + f".tmp{os.getpid()}"
        shutil.rmtree(tmp, ignore_errors=True)
        _split_alevin_by_sample(mtx, row_lbl, n_genes, tmp)
        shutil.rmtree(cache, ignore_errors=True)
        os.rename(tmp, cache)                  # atomic: readers see it complete or not at all
        return cache
    finally:
        try:
            os.unlink(lock)
        except OSError:
            pass


def load_alevin(quant_dir, sample):
    """simpleaf/alevin-fry quant, subset to one sample by row prefix."""
    alv = os.path.join(quant_dir, "alevin")
    mtx = _find(alv, "quants_mat.mtx.gz", "quants_mat.mtx")
    rows = _find(alv, "quants_mat_rows.txt")
    cols = _find(alv, "quants_mat_cols.txt")
    if not (mtx and rows and cols):
        return None
    row_lbl = _open(rows).read().splitlines()
    genes = [c.strip() for c in _open(cols).read().splitlines()]
    pre = sample + "_"
    keep = [i for i, r in enumerate(row_lbl) if r.startswith(pre)]
    if not keep:
        return None
    barcodes = [norm_bc(row_lbl[i]) for i in keep]

    if os.path.getsize(mtx) > SPLIT_ABOVE_BYTES:
        cache = _ensure_split(mtx, alv, row_lbl, len(genes))
        npz = os.path.join(cache, sample + ".npz")
        if not os.path.exists(npz):
            return None
        return sp.load_npz(npz).tocsr(), barcodes, genes

    with _open(mtx) as fh:
        M = scipy.io.mmread(fh).tocsr()            # alevin layout: cells × genes
    return M[keep, :].tocsr(), barcodes, genes


def load_cyto(cyto_root, sample):
    d = os.path.join(cyto_root, "counts", sample)
    return load_10x_dir(d) if os.path.isdir(d) else None


# ---------------------------------------------------------------- gene axis

def probe_name_to_id(probe_csv):
    name2id = {}
    with open(probe_csv) as f:
        rdr = (l for l in f if not l.startswith("#"))
        header = next(rdr).rstrip("\n").split(",")
        gi, ni = header.index("gene_id"), header.index("gene_name")
        for line in rdr:
            p = line.rstrip("\n").split(",")
            if len(p) > max(gi, ni):
                name2id.setdefault(p[ni], p[gi])
    return name2id


def to_gene_id(genes, name2id):
    return [g if g.startswith("ENS") else name2id.get(g, g) for g in genes]


def align(A, a_bc, a_g, B, b_bc, b_g):
    """Restrict both matrices to common barcodes ∩ common genes, same order."""
    bi = {b: i for i, b in enumerate(b_bc)}
    ai = {b: i for i, b in enumerate(a_bc)}
    common_bc = [b for b in a_bc if b in bi]
    agi = {g: i for i, g in enumerate(a_g)}
    bgi = {g: i for i, g in enumerate(b_g)}
    common_g = [g for g in a_g if g in bgi]
    if not common_bc or not common_g:
        return None
    Aa = A[[ai[b] for b in common_bc], :][:, [agi[g] for g in common_g]]
    Bb = B[[bi[b] for b in common_bc], :][:, [bgi[g] for g in common_g]]
    return Aa.tocsr(), Bb.tocsr(), common_bc, common_g


# ---------------------------------------------------------------- layout

def _bench():
    """Layout root (lib/conc_common.sh build_bench_root exports FLEX_LAYOUT_ROOT); KeyError when unset."""
    return os.environ["FLEX_LAYOUT_ROOT"]


DATASETS_TSV = os.environ.get(   # default: the resolved copy config.env writes (config.env DATASETS_TSV)
    "FLEX_DATASETS_TSV",
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "datasets_resolved.tsv"))


def dataset_row(dataset, tsv=None):
    cols = ("id organism plex chemistry geometry cell_bc_whitelist sample_bc_list "
            "probe_csv piscem_index cr_multi_config cr_reference fastq_dir tools").split()
    with open(tsv or DATASETS_TSV) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if p and p[0] == dataset:
                return dict(zip(cols, p))
    raise KeyError(f"{dataset} not in {tsv or DATASETS_TSV}")


def samples_of(dataset, tsv=None):
    """-> [(sample_name, probe_suffix)] from the sample-barcode list actually used.

    Prefers ``sample_bc_list_perprobe.tsv`` when present (see module docstring).
    """
    row = dataset_row(dataset, tsv)
    sbc = row["sample_bc_list"]
    per = os.path.join(os.path.dirname(sbc), "sample_bc_list_perprobe.tsv")
    if os.path.exists(per):
        sbc = per
    out, seen = [], set()
    for line in open(sbc):
        p = line.rstrip("\n").split("\t")
        if len(p) < 3:
            continue
        name, suffix = p[2], p[1]
        if name not in seen:
            seen.add(name)
            out.append((name, suffix))
        elif (name, suffix) not in [(n, s) for n, s in out]:
            # pooled sample spanning several probe barcodes: no unique suffix
            out = [(n, None) if n == name else (n, s) for n, s in out]
    return out


def cellranger_path(dataset, sample):
    """Locate the Cell Ranger filtered matrix for one sample.

    Downloaded 10x per-sample files are ``<prefix>_<sample>_[count_]sample_filtered_
    feature_bc_matrix.h5``; the locally-run mouse dataset writes ``<sample>/`` dirs.
    """
    d = os.path.join(_bench(), "concordance", "cellranger", dataset)
    hit = _find(d, f"*_{sample}_*sample_filtered_feature_bc_matrix.h5")
    if hit:
        return hit
    sub = os.path.join(d, sample)
    if os.path.isdir(sub):
        return sub
    h5s = sorted(glob.glob(os.path.join(d, "*sample_filtered_feature_bc_matrix.h5")))
    return h5s[0] if len(h5s) == 1 else None      # pooled file, split by probe suffix


MIN_SAMPLE_BC_MATCH = 0.5


def sample_bc_match_fraction(quant_dir):
    """matched_reads / mapped_reads for a multiplexed quant, or None if not recorded.

    A wrong sample-barcode list does not fail loudly — it silently routes almost nothing
    to a sample and yields a near-empty matrix (see the v1 sample-
    barcode-list-dropped-7-of-8-variants.md, where 0.36 % of mapped reads matched). Both
    numbers are recorded by simpleaf, so the ratio is checked before a quant is used.
    """
    try:
        with open(os.path.join(quant_dir, "sample_info.json")) as f:
            matched = json.load(f).get("matched_reads")
        with open(os.path.join(quant_dir, "simpleaf_map_info.json")) as f:
            mapped = json.load(f).get("num_mapped")
    except (OSError, ValueError):
        return None
    if not matched or not mapped:
        return None
    return matched / mapped


def simpleaf_quant_dir(dataset, min_match=MIN_SAMPLE_BC_MATCH):
    for c in (os.path.join(_bench(), "concordance", "quant", dataset, "simpleaf", "af_quant"),
              *(os.path.join(_bench(), "out", dataset, f"simpleaf_rep{r}", "af_quant")
                for r in (1, 2, 3))):
        if not _find(os.path.join(c, "alevin"), "quants_mat.mtx*"):
            continue
        frac = sample_bc_match_fraction(c)
        if frac is not None and frac < min_match:
            sys.stderr.write(
                f"[flexio] skipping {c}: only {100 * frac:.2f}% of mapped reads matched a "
                f"sample barcode (threshold {100 * min_match:.0f}%)\n")
            continue
        return c
    return None


def cyto_root(dataset):
    for c in (os.path.join(_bench(), "concordance", "quant", dataset, "cyto"),
              *(os.path.join(_bench(), "out", dataset, f"cyto_rep{r}") for r in (1, 2, 3))):
        if os.path.isdir(os.path.join(c, "counts")):
            return c
    return None

def has_tool(dataset, sample, tool):
    """Cheap existence test for a (dataset, sample, tool) matrix.

    The driver asks this once per sample per tool per round. Answering it by *loading* the
    matrix cost minutes per call on 320k_GEMX (a 15 GB alevin matrix and per-sample cyto
    triplets), so a scheduling round spent half an hour deciding what to run. These checks
    touch metadata only.
    """
    if tool == "cellranger":
        return bool(cellranger_path(dataset, sample))
    if tool == "simpleaf":
        d = simpleaf_quant_dir(dataset)
        if not d:
            return False
        alv = os.path.join(d, "alevin")
        if os.path.exists(os.path.join(alv, "_by_sample", sample + ".npz")):
            return True
        rows = _find(alv, "quants_mat_rows.txt")
        if not rows:
            return False
        pre = (sample + "_").encode()
        with open(rows, "rb") as fh:              # first matching line is enough
            for line in fh:
                if line.startswith(pre):
                    return True
        return False
    if tool == "cyto":
        r = cyto_root(dataset)
        if not r:
            return False
        d = os.path.join(r, "counts", sample)
        return any(os.path.exists(os.path.join(d, n))
                   for n in ("matrix.mtx.gz", "matrix.mtx"))
    raise ValueError(tool)
