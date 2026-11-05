#!/usr/bin/env python3
"""Supp. Figure 1 input: cell-type labels of every cell of one sample, separately for each tool.

The same labelling as the reproduction's step 11 (lib/celltype_annotate.py + lib/cell_labels.py),
applied to each tool's own step-09 per-sample object instead of the simpleaf gene matrix only:
  1  CellTypist Immune_All_Low on the object's raw counts (layer "counts"), normalize_total 1e4 +
     log1p, gene symbols from var.gene_name, majority voting over the object's leiden_res1_0;
  2  the four marker programs of cell_labels.py (erythroid, plasmablast, germinal_centre_B,
     proliferating) scored per leiden_res1_0 cluster; a cluster whose robust z is >= 10 is
     relabelled with the program name (erythroid also needs 5x the median pct_counts_hb).
The step-09 objects carry the common gene axis after its 20-cell gene filter, so the input genes
are a subset of step 11's full simpleaf matrix; all three tools are labelled the same way.

Usage: supp1_tool_labels.py --object H5AD --model PKL --out TSV.GZ
"""
import argparse

import numpy as np
import pandas as pd

PROGRAMS = {   # cell_labels.py PROGRAMS, human symbols
    "erythroid": ["HBA1", "HBA2", "HBB", "ALAS2", "GYPA", "SLC4A1", "CA2"],
    "plasmablast": ["JCHAIN", "XBP1", "PRDM1", "SDC1", "DERL3", "MZB1", "TNFRSF17"],
    "germinal_centre_B": ["AICDA", "GCSAM", "FAS", "MEF2B", "RGS13", "S1PR2"],
    "proliferating": ["MKI67", "TOP2A", "BIRC5", "STMN1", "CCNB1", "UBE2C"],
}
Z_THRESHOLD = 10.0; HB_FOLD_OVER_MEDIAN = 5.0; CLUSTER_KEY = "leiden_res1_0"

ap = argparse.ArgumentParser()
ap.add_argument("--object", required=True); ap.add_argument("--model", required=True)
ap.add_argument("--out", required=True)
o = ap.parse_args()

import anndata as ad
import scanpy as sc
import celltypist

O = ad.read_h5ad(o.object)
A = ad.AnnData(O.layers["counts"].astype(np.float32), obs=O.obs[[CLUSTER_KEY, "pct_counts_hb"]].copy())
A.var_names = O.var["gene_name"].astype(str).values
A.var_names_make_unique()
sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
A.obs[CLUSTER_KEY] = A.obs[CLUSTER_KEY].astype(str)
pred = celltypist.annotate(A, model=o.model, majority_voting=True, over_clustering=CLUSTER_KEY)
ct = pred.predicted_labels.copy()

scores = {}
for p, genes in PROGRAMS.items():
    genes = [g for g in genes if g in A.var_names]
    if genes:
        sc.tl.score_genes(A, genes, score_name=f"_s_{p}", ctrl_size=50)
        scores[p] = A.obs[f"_s_{p}"].values
df = pd.DataFrame(scores, index=A.obs_names); df["cluster"] = A.obs[CLUSTER_KEY].values
per = df.groupby("cluster").mean()
z = per.copy()
for p in per.columns:
    v = per[p].values; mad = np.median(np.abs(v - np.median(v))) * 1.4826
    z[p] = (v - np.median(v)) / (mad if mad > 0 else np.inf)
hb = A.obs.groupby(CLUSTER_KEY)["pct_counts_hb"].mean(); hb_med = float(np.median(hb.values))
remap = {}
for c in z.index:
    called = [p for p in z.columns if z.loc[c, p] >= Z_THRESHOLD]
    if "erythroid" in called and hb.get(c, 0.0) < max(hb_med * HB_FOLD_OVER_MEDIAN, 1e-9):
        called.remove("erythroid")
    remap[c] = "+".join(sorted(called)) if called else ""
ct["cluster"] = A.obs[CLUSTER_KEY].values
ct["marker_program"] = ct["cluster"].map(remap)
ct["label"] = np.where(ct["marker_program"] != "", ct["marker_program"], ct["majority_voting"].astype(str))
ct.index.name = "barcode"
ct[["predicted_labels", "majority_voting", "cluster", "marker_program", "label"]].to_csv(o.out, sep="\t")
print(f"[supp1] {o.object}: {len(ct)} cells, {ct.label.nunique()} labels", flush=True)
