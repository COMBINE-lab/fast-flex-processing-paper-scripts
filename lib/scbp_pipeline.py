#!/usr/bin/env python3
"""
The single-cell best-practices pipeline (https://www.sc-best-practices.org), run
unchanged on each tool's count matrix for one Flex sample, so the resulting clusterings
can be compared between Cell Ranger, simpleaf/alevin-fry and cyto.

Steps, following the book chapter by chapter (parameters are the book's, verbatim):

  preprocessing_visualization/quality_control
    - sc.pp.calculate_qc_metrics(qc_vars=["mt","ribo","hb"], percent_top=[20], log1p=True)
    - MAD outliers: 5 MADs on log1p_total_counts, log1p_n_genes_by_counts,
      pct_counts_in_top_20_genes; 3 MADs on pct_counts_mt OR pct_counts_mt > 8
    - ambient RNA correction with SoupX (soupx_groups from a quick Leiden)
    - sc.pp.filter_genes(min_cells=20)
    - doublets with scDblFinder — ANNOTATED, NOT REMOVED ("We advise to leave identified
      doublets in the dataset for now and inspect doublets during visualization.")
  preprocessing_visualization/normalization
    - scran pooling size factors (computeSumFactors, min.mean=0.1) -> layer
      "scran_normalization"; shifted logarithm kept alongside as "log1p_norm"
  preprocessing_visualization/feature_selection
    - scry::devianceFeatureSelection on raw counts, top 4000 -> var["highly_deviant"]
  preprocessing_visualization/dimensionality_reduction
    - adata.X = layers["scran_normalization"]; highly_variable = highly_deviant
    - sc.pp.pca(svd_solver="arpack", mask_var="highly_variable")
  cellular_structure/clustering
    - sc.pp.neighbors(n_pcs=30)
    - sc.tl.leiden(flavor="igraph", n_iterations=2) at resolutions 0.25, 0.5, 1.0

Two choices the book cannot make for us, because it analyses one matrix and we compare
three:
  * cells — the book starts from a *filtered* (called-cell) matrix. We use Cell Ranger's
    per-sample filtered barcodes for every tool, so differences in the clusterings come
    from quantification and not from three different cell callers. SoupX still gets each
    tool's own raw droplet table for the ambient profile.
  * genes — with --gene-axis common (default) the pipeline runs on the gene_ids shared by
    all available tools, so feature selection sees the same universe; --gene-axis own
    keeps each tool's full gene set.

Two options, both off by default: --remove-doublets (drop scDblFinder doublets before
normalization) and --n-deviant N (deviance feature count, default 4000).

Usage:
  scbp_pipeline.py --dataset 40k_DTC --sample Kidney_Cancer1_BC1_AB1 --tool simpleaf \
      --outdir DIR [--gene-axis common|own] [--threads 8]
"""
import argparse
import json
import logging
import os
import sys

# Cap every numerical backend BEFORE numpy/scanpy/rpy2 load: OpenBLAS, MKL, numexpr and
# numba each default to one thread per core, so a nominally 4-thread run was taking 5-8
# cores (250-380 OS threads) and several in parallel blew the host budget. These variables
# are only read at import time, so they have to be set here rather than in main().
_t = "4"
if "--threads" in sys.argv:
    _t = sys.argv[sys.argv.index("--threads") + 1]
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS",
           "RAYON_NUM_THREADS"):
    os.environ[_v] = _t

# Point rpy2 at THIS environment's R before it is imported: the conda env ships the
# book's Bioconductor packages, while a bare `import rpy2` would find the system R
# (which has none of them), and the user's personal R library would shadow the env's
# with binaries built against a different R.
if "R_HOME" not in os.environ:
    _r_home = os.path.join(sys.prefix, "lib", "R")
    if os.path.isdir(_r_home):
        os.environ["R_HOME"] = _r_home
        os.environ["R_LIBS"] = os.path.join(_r_home, "library")
        os.environ["R_LIBS_USER"] = os.path.join(_r_home, "library")

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix, issparse
from scipy.stats import median_abs_deviation

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flexio

TOOLS = ("cellranger", "simpleaf", "cyto")
SEED = 0          # set from --seed; a second seed on the SAME matrix measures how much of
                  # a tool-vs-tool disagreement is just Leiden/scDblFinder stochasticity


# --------------------------------------------------------------- R bridge

def r_setup(threads):
    # Matrices cross the bridge as dgCMatrix through _to_r_sparse/_from_r_sparse rather
    # than via anndata2ri's AnnData<->SingleCellExperiment conversion, so nothing here
    # depends on that package's (version-dependent) global converter.
    import rpy2.rinterface_lib.callbacks as rcb
    import rpy2.robjects as ro
    rcb.logger.setLevel(logging.ERROR)
    ro.r("suppressPackageStartupMessages({library(Matrix); library(SoupX); "
         "library(scDblFinder); library(scran); library(scry); "
         "library(SingleCellExperiment); library(BiocParallel)})")
    ro.r(f"set.seed({SEED}); register(MulticoreParam({int(threads)}, RNGseed={SEED}))")
    return ro


# --------------------------------------------------------------- book steps

def is_outlier(adata, metric: str, nmads: int):
    """Book verbatim (preprocessing_visualization/quality_control)."""
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def qc_metrics(adata, organism):
    mt_pref = "MT-" if organism == "human" else "mt-"
    rb_pat = ("^RP[SL]", "^Rp[sl]")[organism != "human"]
    hb_pat = ("^HB[^(P)]", "^Hb[^(p)]")[organism != "human"]
    names = adata.var["gene_name"].astype(str)
    adata.var["mt"] = names.str.startswith(mt_pref)
    adata.var["ribo"] = names.str.contains(rb_pat, regex=True)
    adata.var["hb"] = names.str.contains(hb_pat, regex=True)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )
    return adata


def qc_filter(adata):
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
        adata.obs["pct_counts_mt"] > 8
    )
    kept = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    return kept


def soupx_groups(adata):
    adata_pp = adata.copy()
    sc.pp.normalize_total(adata_pp, target_sum=1e4)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="soupx_groups", flavor="igraph",
                 n_iterations=2, directed=False, random_state=SEED)
    return adata_pp.obs["soupx_groups"]


def soupx_correct(ro, adata, adata_raw):
    """SoupX ambient-RNA correction; adata_raw is the tool's table of droplets."""
    import rpy2.robjects as robjects

    groups = soupx_groups(adata)
    adata.obs["soupx_groups"] = groups
    data = adata.X.T.tocsc() if issparse(adata.X) else csr_matrix(adata.X).T.tocsc()
    tod = adata_raw.X.T.tocsc() if issparse(adata_raw.X) else csr_matrix(adata_raw.X).T.tocsc()

    ro.globalenv["data"] = _to_r_sparse(ro, data)
    ro.globalenv["data_tod"] = _to_r_sparse(ro, tod)
    ro.globalenv["genes"] = robjects.StrVector(list(adata.var_names))
    ro.globalenv["cells"] = robjects.StrVector(list(adata.obs_names))
    ro.globalenv["soupx_groups"] = robjects.StrVector([str(g) for g in groups])
    ro.r('''
        rownames(data) = genes; colnames(data) = cells
        rownames(data_tod) = genes
        sc = SoupChannel(data_tod, data, calcSoupProfile = FALSE)
        soupProf = data.frame(row.names = rownames(data),
                              est = rowSums(data)/sum(data), counts = rowSums(data))
        sc = setSoupProfile(sc, soupProf)
        names(soupx_groups) = cells
        sc = setClusters(sc, soupx_groups)
        sc = autoEstCont(sc, doPlot = FALSE)
        out = adjustCounts(sc, roundToInt = TRUE)
        rho = mean(sc$metaData$rho)
    ''')
    out = _from_r_sparse(ro, "out")            # genes × cells
    rho = float(np.asarray(ro.r("rho"))[0])
    adata.layers["counts"] = adata.X.copy()
    adata.layers["soupX_counts"] = csr_matrix(out.T)
    adata.X = adata.layers["soupX_counts"]
    adata.uns["soupx_rho"] = rho
    return adata


def scdblfinder(ro, adata):
    ro.globalenv["data_mat"] = _to_r_sparse(
        ro, adata.X.T.tocsc() if issparse(adata.X) else csr_matrix(adata.X).T.tocsc())
    ro.r(f'''
        set.seed({SEED})
        sce = scDblFinder(SingleCellExperiment(list(counts = data_mat)))
        doublet_score = sce$scDblFinder.score
        doublet_class = as.character(sce$scDblFinder.class)
    ''')
    adata.obs["scDblFinder_score"] = np.asarray(ro.r("doublet_score"))
    adata.obs["scDblFinder_class"] = list(ro.r("doublet_class"))
    return adata


def normalize(ro, adata):
    # shifted logarithm (kept alongside; the book's default recommendation)
    scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
    adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

    # scran pooling size factors — the normalization the clustering chapter's input used
    adata_pp = adata.copy()
    sc.pp.normalize_total(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp, n_comps=15)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="groups", resolution=0.5, flavor="igraph",
                 n_iterations=2, directed=False, random_state=SEED)
    import rpy2.robjects as robjects
    ro.globalenv["data_mat"] = _to_r_sparse(
        ro, adata.X.T.tocsc() if issparse(adata.X) else csr_matrix(adata.X).T.tocsc())
    ro.globalenv["input_groups"] = robjects.StrVector(
        [str(g) for g in adata_pp.obs["groups"]])
    ro.r('''
        size_factors = sizeFactors(computeSumFactors(
            SingleCellExperiment(list(counts = data_mat)),
            clusters = input_groups, min.mean = 0.1, BPPARAM = MulticoreParam()))
    ''')
    adata.obs["size_factors"] = np.asarray(ro.r("size_factors"))
    scran = adata.X / adata.obs["size_factors"].values[:, None]
    adata.layers["scran_normalization"] = csr_matrix(np.log1p(scran))
    return adata


def deviance_features(ro, adata, n_top=4000):
    ro.globalenv["data_mat"] = _to_r_sparse(
        ro, adata.X.T.tocsc() if issparse(adata.X) else csr_matrix(adata.X).T.tocsc())
    ro.r('''
        sce = devianceFeatureSelection(SingleCellExperiment(list(counts = data_mat)),
                                       assay = "counts")
        binomial_deviance = rowData(sce)$binomial_deviance
    ''')
    binomial_deviance = np.asarray(ro.r("binomial_deviance")).ravel()
    idx = binomial_deviance.argsort()[-n_top:]
    mask = np.zeros(adata.var_names.shape, dtype=bool)
    mask[idx] = True
    adata.var["highly_deviant"] = mask
    adata.var["binomial_deviance"] = binomial_deviance
    return adata


def reduce_and_cluster(adata, resolutions=(0.25, 0.5, 1.0)):
    adata.X = adata.layers["scran_normalization"]
    adata.var["highly_variable"] = adata.var["highly_deviant"]
    sc.pp.pca(adata, svd_solver="arpack", mask_var="highly_variable", random_state=SEED)
    sc.pp.neighbors(adata, n_pcs=30, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    for r in resolutions:
        key = "leiden_res" + str(r).replace(".", "_")
        sc.tl.leiden(adata, key_added=key, resolution=r, flavor="igraph",
                     n_iterations=2, random_state=SEED)
    return adata


# --------------------------------------------------------------- sparse bridge

def _to_r_sparse(ro, m):
    """scipy csc -> R dgCMatrix without densifying."""
    import rpy2.robjects as robjects
    from rpy2.robjects import numpy2ri
    m = m.tocsc()
    # R's dgCMatrix requires row indices sorted within each column and no duplicate
    # entries. scipy makes no such guarantee — matrices assembled from COO triplets (as the
    # per-sample split does) can carry unsorted indices, and R then rejects the object with
    # "'i' slot is not increasing within columns" only once it reaches validObject.
    m.sum_duplicates()
    m.sort_indices()
    assert m.has_sorted_indices
    with numpy2ri.converter.context():                 # names must be valid R symbols
        ro.globalenv["zz_i"] = m.indices.astype(np.int32)
        ro.globalenv["zz_p"] = m.indptr.astype(np.int32)
        ro.globalenv["zz_x"] = m.data.astype(float)
        ro.globalenv["zz_d"] = np.array(m.shape, dtype=np.int32)
    return ro.r('new("dgCMatrix", i=as.integer(zz_i), p=as.integer(zz_p), '
                'x=as.numeric(zz_x), Dim=as.integer(zz_d))')


def _from_r_sparse(ro, name):
    from scipy.sparse import csc_matrix
    ro.r(f'{name} <- as(as({name}, "CsparseMatrix"), "dgCMatrix")')
    i = np.asarray(ro.r(f"{name}@i")).astype(np.int64)
    p = np.asarray(ro.r(f"{name}@p")).astype(np.int64)
    x = np.asarray(ro.r(f"{name}@x")).astype(float)
    d = np.asarray(ro.r(f"{name}@Dim")).astype(np.int64)
    return csc_matrix((x, i, p), shape=tuple(d))


# --------------------------------------------------------------- inputs

def build_adata(dataset, sample, tool, gene_axis):
    """-> (adata over Cell Ranger's called cells, adata_raw table of droplets)."""
    import anndata as ad

    row = flexio.dataset_row(dataset)
    name2id = flexio.probe_name_to_id(row["probe_csv"])
    id2name = {}
    for nm, gid in name2id.items():
        id2name.setdefault(gid, nm)
    suffix = dict(flexio.samples_of(dataset)).get(sample)

    cr_path = flexio.cellranger_path(dataset, sample)
    pooled = os.path.isfile(cr_path) and sample not in os.path.basename(cr_path)
    cr = flexio.load_cellranger(cr_path, probe_suffix=suffix if pooled else None)
    if cr is None:
        raise SystemExit(f"no Cell Ranger matrix for {dataset}/{sample}")
    cr_M, cr_bc, cr_g = cr
    cr_g = flexio.to_gene_id(cr_g, name2id)

    loaded = {"cellranger": (cr_M, cr_bc, cr_g)}
    sa_dir = flexio.simpleaf_quant_dir(dataset)
    if sa_dir:
        t = flexio.load_alevin(sa_dir, sample)
        if t:
            loaded["simpleaf"] = (t[0], t[1], flexio.to_gene_id(t[2], name2id))
    cy = flexio.cyto_root(dataset)
    if cy:
        t = flexio.load_cyto(cy, sample)
        if t:
            loaded["cyto"] = (t[0], t[1], flexio.to_gene_id(t[2], name2id))
    if tool not in loaded:
        raise SystemExit(f"{tool} matrix unavailable for {dataset}/{sample}")

    if gene_axis == "common":
        genes = set(loaded["cellranger"][2])
        for _, _, g in loaded.values():
            genes &= set(g)
        keep_genes = [g for g in loaded["cellranger"][2] if g in genes]
    else:
        keep_genes = None

    def to_adata(M, bcs, gs, cells):
        gi = {g: i for i, g in enumerate(gs)}
        bi = {b: i for i, b in enumerate(bcs)}
        cols = [gi[g] for g in (keep_genes or gs) if g in gi]
        gene_ids = [g for g in (keep_genes or gs) if g in gi]
        rows = [bi[b] for b in cells if b in bi]
        used = [b for b in cells if b in bi]
        A = ad.AnnData(M[rows, :][:, cols].tocsr().astype(np.float32))
        A.obs_names = used
        A.var_names = gene_ids
        A.var["gene_name"] = [id2name.get(g, g) for g in gene_ids]
        return A

    M, bcs, gs = loaded[tool]
    adata = to_adata(M, bcs, gs, cr_bc)                    # Cell Ranger's called cells
    raw_cells = bcs if tool != "cellranger" else None
    if tool == "cellranger":
        raw_path = _cellranger_raw(dataset, sample)
        if raw_path:
            rM, rbc, rg = flexio.load_cellranger(raw_path,
                                                 probe_suffix=suffix if pooled else None)
            adata_raw = to_adata(rM, rbc, flexio.to_gene_id(rg, name2id), rbc)
        else:
            adata_raw = adata.copy()
    else:
        adata_raw = to_adata(M, bcs, gs, raw_cells)        # every barcode the tool emitted
    return adata, adata_raw


def _cellranger_raw(dataset, sample):
    d = os.path.join(flexio._bench(), "concordance", "cellranger", dataset)
    hit = flexio._find(d, f"*_{sample}_*sample_raw_feature_bc_matrix.h5")
    if hit:
        return hit
    sub = os.path.join(d, sample, "sample_raw_feature_bc_matrix.h5")
    if os.path.exists(sub):
        return sub
    import glob
    h5s = sorted(glob.glob(os.path.join(d, "*sample_raw_feature_bc_matrix.h5")))
    return h5s[0] if len(h5s) == 1 else None


# --------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--tool", required=True, choices=TOOLS)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--gene-axis", default="common", choices=("common", "own"))
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--remove-doublets", action="store_true",
                    help="drop scDblFinder_class == 'doublet' cells right after scDblFinder, before "
                         "normalization / feature selection / PCA / Leiden (default: annotate only, the book's advice)")
    ap.add_argument("--n-deviant", type=int, default=4000,
                    help="number of most deviant genes kept for PCA (default 4000, the book's)")
    ap.add_argument("--seed", type=int, default=0,
                    help="re-run identically with a different seed to get the pipeline's "
                         "own run-to-run variation (the ceiling any tool pair can reach)")
    args = ap.parse_args()

    global SEED
    SEED = args.seed

    sc.settings.n_jobs = args.threads
    sc.settings.verbosity = 1

    row = flexio.dataset_row(args.dataset)
    ro = r_setup(args.threads)

    adata, adata_raw = build_adata(args.dataset, args.sample, args.tool, args.gene_axis)
    n_in = adata.n_obs
    adata = qc_metrics(adata, row["organism"])
    adata = qc_filter(adata)
    adata_raw = adata_raw[:, adata.var_names].copy()
    adata = soupx_correct(ro, adata, adata_raw)
    sc.pp.filter_genes(adata, min_cells=20)
    adata = scdblfinder(ro, adata)
    n_qc = adata.n_obs
    n_doublets = int((adata.obs["scDblFinder_class"] == "doublet").sum())
    if args.remove_doublets:
        adata = adata[adata.obs["scDblFinder_class"] != "doublet"].copy()
    adata = normalize(ro, adata)
    adata = deviance_features(ro, adata, n_top=args.n_deviant)
    adata = reduce_and_cluster(adata)

    tag = f"{args.dataset}__{args.sample}__{args.tool}__{args.gene_axis}"
    if args.seed:
        tag += f"__seed{args.seed}"
    os.makedirs(args.outdir, exist_ok=True)
    keep = [c for c in adata.obs.columns
            if c.startswith("leiden_res") or c in
            ("total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo",
             "pct_counts_hb", "pct_counts_in_top_20_genes", "size_factors",
             "scDblFinder_score", "scDblFinder_class", "soupx_groups")]
    obs = adata.obs[keep].copy()
    obs.index.name = "barcode"
    obs.to_csv(os.path.join(args.outdir, tag + ".obs.tsv.gz"), sep="\t")
    adata.write_h5ad(os.path.join(args.outdir, tag + ".h5ad"), compression="gzip")

    summary = {
        "dataset": args.dataset, "sample": args.sample, "tool": args.tool,
        "gene_axis": args.gene_axis, "seed": args.seed,
        "n_cells_in": int(n_in), "n_cells_after_qc": int(n_qc),
        "n_genes_after_filter": int(adata.n_vars),
        "soupx_rho": float(adata.uns.get("soupx_rho", float("nan"))),
        "n_doublets": n_doublets,
        "n_clusters": {k: int(adata.obs[k].nunique())
                       for k in adata.obs if k.startswith("leiden_res")},
    }
    if args.remove_doublets or args.n_deviant != 4000:      # variant runs only: defaults keep the paper's file
        summary.update({"remove_doublets": bool(args.remove_doublets), "n_doublets_removed":
                        n_doublets if args.remove_doublets else 0, "n_cells_final": int(adata.n_obs),
                        "n_deviant": args.n_deviant})
    with open(os.path.join(args.outdir, tag + ".summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
