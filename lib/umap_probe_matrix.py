#!/usr/bin/env python3
"""Pooled UMAP of one library's donors, for Figure 2 panel b and Supp. Figure 2.

The embedding comes from the gene matrix (probes summed per gene), so the layout is the
ordinary gene-level analysis and the color is what only the probe matrix can provide.
A single probe carries under one UMI per cell, so the panel script pools each cell's probe
values over its nearest neighbors in the same PCA space.

Donors default to every probe matrix of the library. Paths, k and the seed are arguments. The
written h5ad keeps only obs, PCA and UMAP; the panel script re-reads the probe matrices.

Usage:
  umap_probe_matrix.py --probe-root DIR --label-root DIR --dataset DS --out H5AD
                       [--samples S ...] [--n-neighbors 15] [--seed 0]
"""
import argparse
import os
import sys

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from probe_instruments import probe_paths, read_labels, write_doublet_log


def load_labelled(probe_root, label_root, dataset, samples=None, drop_doublets=False):
    """Concatenate the labelled cells of each sample; obs_names become <sample>_<barcode>."""
    import anndata as ad
    parts, var = [], None
    for d, path in probe_paths(probe_root, dataset, samples):
        a = ad.read_h5ad(path)
        lab = read_labels(label_root, dataset, d, drop_doublets)
        keep = [b for b in a.obs_names if b in lab.index]
        if not keep:
            raise SystemExit(f"{d}: no barcode shared between matrix and labels")
        a = a[keep].copy()
        a.obs["label"] = lab.loc[keep].values
        a.obs["donor"] = d
        a.obs_names = [f"{d}_{b}" for b in a.obs_names]
        if var is None:
            var = a.var.copy()
        elif not a.var.index.equals(var.index):
            raise SystemExit(f"{d}: probe axis differs from the first sample")
        parts.append(a)
        print(d, a.shape, flush=True)
    a = ad.concat(parts, join="outer", merge="first")
    a.var = var.loc[a.var_names]
    return a


def main():
    import anndata as ad
    import scanpy as sc
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe-root", required=True); ap.add_argument("--label-root", required=True)
    ap.add_argument("--dataset", required=True); ap.add_argument("--samples", nargs="+")
    ap.add_argument("--drop-doublets", action="store_true")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-neighbors", type=int, default=15)
    ap.add_argument("--seed", type=int, default=0)
    o = ap.parse_args()

    a = load_labelled(o.probe_root, o.label_root, o.dataset, o.samples, o.drop_doublets)
    if o.drop_doublets:
        write_doublet_log(os.path.splitext(o.out)[0] + ".doublets_removed.tsv")
    print("pooled", a.shape, flush=True)
    gid = a.var["gene_id"].astype(str).values
    genes, inv = np.unique(gid, return_inverse=True)
    M = sp.csr_matrix((np.ones(len(inv)), (np.arange(len(inv)), inv)),
                      shape=(len(inv), len(genes)))
    G = ad.AnnData(sp.csr_matrix(a.X @ M), obs=a.obs.copy())
    G.var_names = genes
    G.obs["depth"] = np.asarray(G.X.sum(axis=1)).ravel()

    sc.pp.normalize_total(G, target_sum=1e4); sc.pp.log1p(G)
    sc.pp.highly_variable_genes(G, n_top_genes=2000)
    sc.pp.pca(G, n_comps=50, mask_var="highly_variable", random_state=o.seed)
    sc.pp.neighbors(G, n_neighbors=o.n_neighbors, random_state=o.seed)
    sc.tl.umap(G, random_state=o.seed)
    donors = sorted(set(G.obs["donor"]))
    out = ad.AnnData(obs=G.obs[["label", "donor", "depth"]].copy(),
                     obsm={"X_umap": G.obsm["X_umap"], "X_pca": G.obsm["X_pca"]},
                     uns={"donors": donors, "n_neighbors": o.n_neighbors, "seed": o.seed})
    out.write_h5ad(o.out)
    print(f"wrote {o.out}: {out.n_obs} cells, {len(donors)} donors", flush=True)


if __name__ == "__main__":
    main()
