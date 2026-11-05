#!/usr/bin/env python3
"""
Count-matrix CONCORDANCE, Cell Ranger as reference (benchmark axis 2).

Metric set is the alevin-fry paper's (Nat. Methods 2022, supp. S8 / `notebooks/
starsolo_sim_analysis.ipynb` in COMBINE-lab/alevin-fry-paper-scripts):
  mean & median per-cell Spearman, MARD (drop-NA and NA=0), mean rFP/cell, mean rFN/cell,
plus the cyto-preprint style count concordance % and cell/gene overlap, and per-gene
correlation of summed counts.

The cell set is Cell Ranger's *filtered* (called) barcodes for the sample, applied
unchanged to every other tool — one fixed cell set for all tools, so the numbers measure
quantification differences and not each tool's own cell calling. Genes are matched on
gene_id (cyto's gene-name axis harmonized via the probe-set CSV). Matrices are compared
RAW: no normalization, no gene/cell filtering.

Usage:
  compare_consistency.py --dataset 40k_DTC --sample Kidney_Cancer1_BC1_AB1 --out cmp.json
  compare_consistency.py --dataset ... --sample ... --cellranger PATH --simpleaf DIR --cyto DIR
"""
import argparse
import json
import os
import sys

import numpy as np
import scipy.sparse as sp
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flexio


def metrics(X, Y, spearman_sample=None, seed=0):
    """X=reference, Y=test — both csr cells×genes, aligned, RAW.

    Computed sparsely: the raw matrices are far too large to densify.
    """
    X = X.tocsr().astype(float)
    Y = Y.tocsr().astype(float)
    n_cells, n_genes = X.shape

    # --- per-gene summed-count correlation ---
    xg = np.asarray(X.sum(axis=0)).ravel()
    yg = np.asarray(Y.sum(axis=0)).ravel()
    gene_pearson = float(pearsonr(np.log1p(xg), np.log1p(yg))[0]) if n_genes > 2 else None
    gene_spear = float(spearmanr(xg, yg).correlation) if n_genes > 2 else None

    # --- L1 / count concordance ---
    D = abs(X - Y)
    l1 = float(D.sum())
    tot = float(X.sum() + Y.sum())
    concordance = float(100.0 * (1.0 - l1 / tot)) if tot > 0 else None

    # --- MARD (S8 eq.1 drop-NA, eq.3 NA=0); sparse over the union of nonzeros ---
    Mx = X.maximum(Y).tocsr()
    nnz_union = Mx.nnz
    if nnz_union > 0:
        Mxinv = Mx.copy()
        Mxinv.data = 1.0 / Mxinv.data
        ratio_sum = float(D.multiply(Mxinv).sum())
        mard_dropna = ratio_sum / nnz_union
        mard_na0 = ratio_sum / (n_cells * n_genes)
    else:
        mard_dropna = mard_na0 = None

    # --- rFP/rFN per cell (S8 eq.4,5), relative to genes expressed in the reference ---
    Xb = (X > 0).astype(float)
    Yb = (Y > 0).astype(float)
    both = Xb.multiply(Yb)
    nx = np.asarray(Xb.sum(axis=1)).ravel()
    ny = np.asarray(Yb.sum(axis=1)).ravel()
    nb = np.asarray(both.sum(axis=1)).ravel()
    fp, fn = ny - nb, nx - nb
    valid = nx > 0
    rfp = float((fp[valid] / nx[valid]).mean()) if valid.any() else None
    rfn = float((fn[valid] / nx[valid]).mean()) if valid.any() else None

    # --- per-cell Spearman over all common cells (optionally subsampled) ---
    idx = np.arange(n_cells)
    sampled = spearman_sample is not None and n_cells > spearman_sample
    if sampled:
        idx = np.random.default_rng(seed).choice(n_cells, spearman_sample, replace=False)
    sps = []
    for i in idx:
        xi = X.getrow(i).toarray().ravel()
        yi = Y.getrow(i).toarray().ravel()
        if xi.any() and yi.any():
            r = spearmanr(xi, yi).correlation
            if not np.isnan(r):
                sps.append(r)
    return {
        "n_common_cells": int(n_cells), "n_common_genes": int(n_genes),
        "mean_cell_spearman": float(np.mean(sps)) if sps else None,
        "median_cell_spearman": float(np.median(sps)) if sps else None,
        "cell_spearman_n": len(sps), "cell_spearman_sampled": bool(sampled),
        "MARD_dropNA": mard_dropna, "MARD_NA0": mard_na0,
        "mean_rFP_per_cell": rfp, "mean_rFN_per_cell": rfn,
        "gene_pearson_log1p": gene_pearson, "gene_spearman": gene_spear,
        "count_concordance_pct": concordance,
        "total_umi_ref": float(X.sum()), "total_umi_test": float(Y.sum()),
        "total_umi_ratio_test_over_ref": float(Y.sum() / X.sum()) if X.sum() else None,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--probe-csv", help="default: from datasets_resolved.tsv")
    ap.add_argument("--cellranger", help="filtered .h5 or per_sample_outs dir")
    ap.add_argument("--simpleaf", help="af_quant dir")
    ap.add_argument("--cyto", help="cyto output root (expects counts/<sample>)")
    ap.add_argument("--spearman-sample", type=int, default=None,
                    help="cap the per-cell Spearman estimator at N cells (default: all)")
    args = ap.parse_args()

    row = flexio.dataset_row(args.dataset)
    probe_csv = args.probe_csv or row["probe_csv"]
    name2id = flexio.probe_name_to_id(probe_csv)

    suffix = dict(flexio.samples_of(args.dataset)).get(args.sample)
    cr_path = args.cellranger or flexio.cellranger_path(args.dataset, args.sample)
    if not cr_path:
        sys.exit(f"no Cell Ranger matrix for {args.dataset}/{args.sample}")
    # a per-sample file is already one probe barcode; only split a pooled file
    pooled = os.path.isfile(cr_path) and args.sample not in os.path.basename(cr_path)
    cr = flexio.load_cellranger(cr_path, probe_suffix=suffix if pooled else None)
    if cr is None:
        sys.exit(f"could not load Cell Ranger matrix {cr_path}")
    cr_M, cr_bc, cr_g = cr
    cr_g = flexio.to_gene_id(cr_g, name2id)
    cr_bc_set = set(cr_bc)

    result = {"dataset": args.dataset, "sample": args.sample,
              "cellranger_matrix": cr_path,
              "cell_set": "cellranger_filtered_barcodes (applied to all tools)",
              "n_cellranger_filtered_cells": len(cr_bc),
              "n_cellranger_genes": len(cr_g), "pairs": {}}

    tools = {}
    sa_dir = args.simpleaf or flexio.simpleaf_quant_dir(args.dataset)
    if sa_dir:
        t = flexio.load_alevin(sa_dir, args.sample)
        if t:
            tools["simpleaf"] = (t[0], t[1], flexio.to_gene_id(t[2], name2id))
            result["simpleaf_quant"] = sa_dir
    cy_root = args.cyto or flexio.cyto_root(args.dataset)
    if cy_root:
        t = flexio.load_cyto(cy_root, args.sample)
        if t:
            tools["cyto"] = (t[0], t[1], flexio.to_gene_id(t[2], name2id))
            result["cyto_root"] = cy_root

    for name, (M, bc, g) in tools.items():
        al = flexio.align(cr_M, cr_bc, cr_g, M, bc, g)
        if al is None:
            result["pairs"][f"cellranger_vs_{name}"] = {"error": "no common cells/genes"}
            continue
        X, Y, cbc, cg = al
        m = metrics(X, Y, spearman_sample=args.spearman_sample)
        m["test_total_cells"] = len(bc)
        m["cell_overlap_jaccard"] = len(cbc) / len(cr_bc_set | set(bc))
        m["gene_overlap_jaccard"] = len(cg) / len(set(cr_g) | set(g))
        m["cellranger_cells_recovered_pct"] = 100.0 * len(cbc) / len(cr_bc)
        result["pairs"][f"cellranger_vs_{name}"] = m

    if "simpleaf" in tools and "cyto" in tools:
        def restrict(M, bc, g):
            rows = [i for i, b in enumerate(bc) if b in cr_bc_set]
            return M[rows, :].tocsr(), [bc[i] for i in rows], g
        al = flexio.align(*restrict(*tools["simpleaf"]), *restrict(*tools["cyto"]))
        if al:
            X, Y, _, cg = al
            m = metrics(X, Y, spearman_sample=args.spearman_sample)
            m["gene_overlap_jaccard"] = len(cg) / len(
                set(tools["simpleaf"][2]) | set(tools["cyto"][2]))
            result["pairs"]["simpleaf_vs_cyto"] = m

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    json.dump(result, open(args.out, "w"), indent=2)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
