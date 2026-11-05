#!/usr/bin/env python3
"""
Validate and export a probe-granularity Flex count matrix.

`alevin-fry quant` run with a t2g whose first two columns are both the probe id yields one
feature per probe instead of one per gene. Nothing else about the pipeline changes: the same
collated RAD, the same permit list, the same `cr-like` resolution. So the probe matrix must
aggregate back to the campaign's gene matrix, and `validate` is the gate that proves it did
before any downstream analysis runs.

`export` then writes the called-cell matrices the analysis actually uses: one h5ad per
(dataset, sample) over Cell Ranger's filtered barcodes — the same cell definition
`scbp_pipeline.py` uses for every tool — with each probe's re-derived annotation
(`probe_annot/<panel>.tsv`) carried in `.var`.

Both passes stream: 320k_GEMX's probe matrix holds ~10^9 entries, far more than fits in RAM.

The validate gate is STRUCTURAL only
(|loss| < --loss-max, identical non-zero gene set, log-corr > --logcorr-min). The per-cell check
(`cells_losing_over_1pct == 0`) is recorded as `PASS_strict` but does not gate: at --umi-edit-dist 1
the Hamming-1 collapse is per feature, so two probes of one gene keep UMIs the gene-level collapse
merges (probe sum >= gene sum).
`flexio` is imported from this directory (lib/flexio.py).

Usage:
  probe_matrix.py validate --probe-quant DIR --gene-quant DIR [--json OUT] [--loss-max F] [--logcorr-min F]
  probe_matrix.py export   --probe-quant DIR --dataset DS --annot TSV --outdir DIR --cr-dir DIR
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flexio

CHUNK = 20_000_000


def _dims(mtx):
    with open(mtx) as f:
        for line in f:
            if not line.startswith("%"):
                return tuple(int(x) for x in line.split())
    raise ValueError(f"no dimension line in {mtx}")


def _margins(mtx, n_rows, n_cols):
    """Streamed (row_sums, col_sums, total) of a MatrixMarket file."""
    rs = np.zeros(n_rows + 1)
    cs = np.zeros(n_cols + 1)
    with open(mtx) as f:
        for line in f:
            if not line.startswith("%"):
                break                                    # the dims line, consumed
        for chunk in pd.read_csv(f, sep=r"\s+", header=None, chunksize=CHUNK,
                                 names=("r", "c", "v"),
                                 dtype={"r": np.int64, "c": np.int64, "v": np.float64}):
            np.add.at(rs, chunk.r.to_numpy(), chunk.v.to_numpy())
            np.add.at(cs, chunk.c.to_numpy(), chunk.v.to_numpy())
    return rs[1:], cs[1:], rs.sum()


def validate(probe_quant, gene_quant, loss_max=0.005, logcorr_min=0.999):
    pa = os.path.join(probe_quant, "alevin")
    ga = os.path.join(gene_quant, "alevin")
    p_rows = open(f"{pa}/quants_mat_rows.txt").read().split()
    p_cols = open(f"{pa}/quants_mat_cols.txt").read().split()
    g_rows = open(f"{ga}/quants_mat_rows.txt").read().split()
    g_cols = open(f"{ga}/quants_mat_cols.txt").read().split()

    # Check the cell axes are the same thing BEFORE streaming ~10^9 entries twice. A quant run
    # against a directory lacking collation_manifest.bin emits bare 16 nt labels instead of
    # <sample>_<barcode>; the matrix looks well-formed and joins to nothing, so it has to be
    # caught by name rather than discovered five minutes later as a KeyError.
    if len(p_rows) != len(g_rows):
        raise SystemExit(f"cell axes differ in length: probe {len(p_rows)} vs gene "
                         f"{len(g_rows)} -- these are not the same quantification")
    shared = len(set(p_rows) & set(g_rows))
    if shared < 0.99 * len(g_rows):
        raise SystemExit(
            f"cell axes do not match: only {shared} of {len(g_rows)} row labels shared.\n"
            f"  probe example: {p_rows[0]!r}\n  gene  example: {g_rows[0]!r}\n"
            f"Bare barcodes on one side mean the quant directory had no "
            f"collation_manifest.bin and the sample assignment was lost.")

    pr, pc, ptot = _margins(f"{pa}/quants_mat.mtx", *_dims(f"{pa}/quants_mat.mtx")[:2])
    gr, gc, gtot = _margins(f"{ga}/quants_mat.mtx", *_dims(f"{ga}/quants_mat.mtx")[:2])

    # features: probe column sums summed per gene, against the gene matrix
    gene_of = np.array([p.split("|")[0] for p in p_cols])
    agg = pd.Series(pc).groupby(gene_of).sum()
    j = pd.concat([agg.rename("probe"), pd.Series(gc, index=g_cols).rename("gene")],
                  axis=1).fillna(0.0)
    rel = np.divide(j.probe - j.gene, j.gene,
                    out=np.zeros(len(j)), where=j.gene.to_numpy() > 0)
    hi = j.gene.to_numpy() >= 1000

    # cells: aligned by barcode label, which is identical in both matrices
    gi = {b: i for i, b in enumerate(g_rows)}
    order = np.array([gi[b] for b in p_rows])
    cell_g = gr[order]
    cell_rel = np.divide(pr - cell_g, cell_g,
                         out=np.zeros(len(pr)), where=cell_g > 0)
    ok_cell = cell_g >= 200

    m = {
        "n_probes": len(p_cols),
        "n_genes": len(g_cols),
        "n_barcodes": len(p_rows),
        "probe_total_umi": float(ptot),
        "gene_total_umi": float(gtot),
        "total_umi_loss_frac": float((gtot - ptot) / gtot) if gtot else None,
        "genes_nonzero_probe": int((j.probe > 0).sum()),
        "genes_nonzero_gene": int((j.gene > 0).sum()),
        "gene_rel_diff_median": float(np.median(np.abs(rel))),
        "gene_rel_diff_p99": float(np.quantile(np.abs(rel), 0.99)),
        "gene_rel_diff_max": float(np.abs(rel).max()),
        "gene_rel_diff_max_ge1000umi": float(np.abs(rel[hi]).max()) if hi.any() else None,
        "gene_logcorr": float(np.corrcoef(np.log1p(j.probe[j.gene > 0]),
                                          np.log1p(j.gene[j.gene > 0]))[0, 1]),
        "cells_ge200umi": int(ok_cell.sum()),
        "cell_rel_diff_median": float(np.median(np.abs(cell_rel[ok_cell]))),
        "cell_rel_diff_max": float(np.abs(cell_rel[ok_cell]).max()),
        "cells_losing_over_1pct": int((np.abs(cell_rel[ok_cell]) > 0.01).sum()),
    }
    # The gate: probe granularity must be a decomposition, not a different measurement.
    m["PASS"] = bool(
        abs(m["total_umi_loss_frac"]) < loss_max
        and m["genes_nonzero_probe"] == m["genes_nonzero_gene"]
        and m["gene_logcorr"] > logcorr_min)
    m["PASS_strict"] = bool(m["PASS"] and m["cells_losing_over_1pct"] == 0)
    m["gate"] = {"loss_max": loss_max, "logcorr_min": logcorr_min, "per_cell_check_gates": False}
    return m


def _called_cells(dataset, cr_dir):
    """{sample: [16 nt barcodes]} -- Cell Ranger's called cells, the analysis cell definition.

    Reproduction-package change: read from our Cell Ranger 10.1.0 run,
    <cr_dir>/<sample>/sample_filtered_feature_bc_matrix.h5 (one file per sample), and fail when a
    sample's file is missing instead of skipping it."""
    out = {}
    for sample, _ in flexio.samples_of(dataset):
        h5 = os.path.join(cr_dir, sample, "sample_filtered_feature_bc_matrix.h5")
        if not os.path.isfile(h5):
            raise SystemExit(f"{dataset}: no Cell Ranger filtered matrix {h5}")
        out[sample] = flexio.load_cellranger(h5)[1]
    return out


def _stream_called(probe_quant, called):
    """One pass over quants_mat.mtx keeping only called-cell rows -> {sample: (csr, bcs)}.

    `flexio.load_alevin` builds a per-sample cache covering *every* barcode the quant emitted,
    which for a 16-plex library is 2.6 M barcodes and ~6e8 entries -- tens of GB of scratch and
    many minutes, when the analysis only ever uses the ~130 k called cells. This keeps the ~5 %
    of entries that belong to a called cell and nothing else.
    """
    alv = os.path.join(probe_quant, "alevin")
    row_lbl = open(os.path.join(alv, "quants_mat_rows.txt")).read().split()
    n_cols = len(open(os.path.join(alv, "quants_mat_cols.txt")).read().split())

    want = {}                                     # row index -> (sample, local row)
    order = {}
    for sample, bcs in called.items():
        order[sample] = []
    pos = {}
    for gi, lbl in enumerate(row_lbl):
        smp = lbl.rsplit("_", 1)[0]
        if smp in called:
            pos.setdefault(smp, {})[flexio.norm_bc(lbl)] = gi
    for sample, bcs in called.items():
        m = pos.get(sample, {})
        for b in bcs:
            gi = m.get(b)
            if gi is not None:
                want[gi] = (sample, len(order[sample]))
                order[sample].append(b)

    # Sample identity is carried as an integer code, not a string: with ~6e8 entries an object
    # array turns every chunk's grouping into Python-object comparisons and dominates runtime.
    names = list(order)
    code = {s: i for i, s in enumerate(names)}
    keep = np.zeros(len(row_lbl), bool)
    keep[list(want)] = True
    smp_of = np.full(len(row_lbl), -1, dtype=np.int32)
    loc_of = np.full(len(row_lbl), -1, dtype=np.int64)
    for gi, (s, lr) in want.items():
        smp_of[gi] = code[s]
        loc_of[gi] = lr

    acc = {s: ([], [], []) for s in names}
    with open(os.path.join(alv, "quants_mat.mtx")) as f:
        for line in f:
            if not line.startswith("%"):
                break
        for ch in pd.read_csv(f, sep=r"\s+", header=None, chunksize=CHUNK,
                              names=("r", "c", "v"), engine="c",
                              dtype={"r": np.int64, "c": np.int32, "v": np.float32}):
            gi = ch.r.to_numpy() - 1
            m = keep[gi]
            if not m.any():
                continue
            gi = gi[m]
            cc = ch.c.to_numpy()[m] - 1
            vv = ch.v.to_numpy()[m]
            ss = smp_of[gi]
            ll = loc_of[gi]
            srt = np.argsort(ss, kind="stable")           # one pass, contiguous per sample
            ss, ll, cc, vv = ss[srt], ll[srt], cc[srt], vv[srt]
            bounds = np.searchsorted(ss, np.arange(len(names) + 1))
            for i, s in enumerate(names):
                a, b = bounds[i], bounds[i + 1]
                if b > a:
                    acc[s][0].append(ll[a:b]); acc[s][1].append(cc[a:b]); acc[s][2].append(vv[a:b])

    out = {}
    for s, (rr, cc, vv) in acc.items():
        if not rr:
            continue
        M = sp.coo_matrix((np.concatenate(vv),
                           (np.concatenate(rr), np.concatenate(cc))),
                          shape=(len(order[s]), n_cols)).tocsr()
        out[s] = (M, order[s])
    return out


def export(probe_quant, dataset, annot_tsv, outdir, cr_dir):
    import anndata as ad

    ann = pd.read_csv(annot_tsv, sep="\t").set_index("probe_id")
    os.makedirs(outdir, exist_ok=True)
    probes = open(os.path.join(probe_quant, "alevin",
                               "quants_mat_cols.txt")).read().split()
    called = _called_cells(dataset, cr_dir)
    print(f"[export] {len(called)} samples, "
          f"{sum(len(v) for v in called.values())} called cells", flush=True)
    streamed = _stream_called(probe_quant, called)

    written = {}
    for sample in called:
        if sample not in streamed:
            print(f"[export] {sample}: no probe matrix rows, skipped", flush=True)
            continue
        M, used = streamed[sample]
        A = ad.AnnData(M.astype(np.float32))
        A.obs_names = used
        A.obs["sample"] = sample
        A.var_names = probes
        for c in ("gene_id", "gene_name", "region_csv", "genomic_class",
                  "n_tx_hit", "n_tx_gene", "constitutive", "tx_set"):
            if c in ann.columns:
                A.var[c] = ann[c].reindex(probes).values
        out = os.path.join(outdir, f"{sample}.probe.h5ad")
        A.write_h5ad(out, compression="gzip")
        written[sample] = {"path": out, "n_cells": A.n_obs, "n_probes": A.n_vars,
                           "n_called": len(called[sample]), "nnz": int(A.X.nnz)}
        print(f"[export] {sample}: {A.n_obs}/{len(called[sample])} called cells x "
              f"{A.n_vars} probes, nnz={A.X.nnz}", flush=True)
    return written


def _export_legacy(probe_quant, dataset, annot_tsv, outdir):
    """Original per-sample loader via flexio's cache; kept as a cross-check path."""
    import anndata as ad

    ann = pd.read_csv(annot_tsv, sep="\t").set_index("probe_id")
    os.makedirs(outdir, exist_ok=True)
    written = {}
    for sample, suffix in flexio.samples_of(dataset):
        cr_path = flexio.cellranger_path(dataset, sample)
        if not cr_path:
            print(f"[export] {sample}: no Cell Ranger matrix, skipped", flush=True)
            continue
        pooled = os.path.isfile(cr_path) and sample not in os.path.basename(cr_path)
        cr = flexio.load_cellranger(cr_path, probe_suffix=suffix if pooled else None)
        if cr is None:
            print(f"[export] {sample}: Cell Ranger load failed, skipped", flush=True)
            continue
        called = cr[1]                                   # 16 nt cores, the cell definition

        loaded = flexio.load_alevin(probe_quant, sample)
        if loaded is None:
            print(f"[export] {sample}: no probe matrix rows, skipped", flush=True)
            continue
        M, bcs, probes = loaded
        bi = {b: i for i, b in enumerate(bcs)}
        rows = [bi[b] for b in called if b in bi]
        used = [b for b in called if b in bi]

        A = ad.AnnData(M[rows, :].tocsr().astype(np.float32))
        A.obs_names = used
        A.obs["sample"] = sample
        A.var_names = probes
        for c in ("gene_id", "gene_name", "region_csv", "genomic_class",
                  "n_tx_hit", "n_tx_gene", "constitutive", "tx_set"):
            if c in ann.columns:
                A.var[c] = ann[c].reindex(probes).values
        out = os.path.join(outdir, f"{sample}.probe.h5ad")
        A.write_h5ad(out, compression="gzip")
        written[sample] = {"path": out, "n_cells": A.n_obs, "n_probes": A.n_vars,
                           "n_called": len(called), "nnz": int(A.X.nnz)}
        print(f"[export] {sample}: {A.n_obs}/{len(called)} called cells x {A.n_vars} "
              f"probes, nnz={A.X.nnz}", flush=True)
    return written


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate")
    v.add_argument("--probe-quant", required=True)
    v.add_argument("--gene-quant", required=True)
    v.add_argument("--json")
    v.add_argument("--loss-max", type=float, default=0.005)
    v.add_argument("--logcorr-min", type=float, default=0.999)
    e = sub.add_parser("export")
    e.add_argument("--probe-quant", required=True)
    e.add_argument("--dataset", required=True)
    e.add_argument("--annot", required=True)
    e.add_argument("--outdir", required=True)
    e.add_argument("--cr-dir", required=True, help="<CR per_sample_outs>/<sample>/sample_filtered_feature_bc_matrix.h5")
    e.add_argument("--json")
    a = ap.parse_args()

    if a.cmd == "validate":
        m = validate(a.probe_quant, a.gene_quant, a.loss_max, a.logcorr_min)
        for k, val in m.items():
            print(f"{k:32s} {val}")
        if a.json:
            json.dump(m, open(a.json, "w"), indent=2)
        sys.exit(0 if m["PASS"] else 1)

    res = export(a.probe_quant, a.dataset, a.annot, a.outdir, a.cr_dir)
    if a.json:
        json.dump(res, open(a.json, "w"), indent=2)
    sys.exit(0 if res else 1)


if __name__ == "__main__":
    main()
