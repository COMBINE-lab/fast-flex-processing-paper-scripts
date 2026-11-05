#!/usr/bin/env python3
"""
Final per-cell label table: CellTypist, refined by marker programs it cannot represent.

`Immune_All_Low` is the only CellTypist model covering these tissues, and it is immune-only:
it has no erythroid category at all, so mouse spleen erythroblasts -- one of the few genuinely
differentiating populations in this benchmark -- are silently absorbed into whatever immune
label is nearest. Majority voting also washes out the low-frequency states that matter here
(`Proliferative germinal center B cells` survives in the raw labels but not the vote).

So four marker programs are scored per Leiden cluster and a cluster is relabelled when its
score stands clearly above the other clusters'. The programs are exactly the ones the model
cannot reach or reliably vote through:

  erythroid      the model has no such class
  plasmablast    terminal B differentiation, low frequency
  germinal_centre_B  active affinity maturation, low frequency
  proliferating  cycling cells of any lineage

Both layers are written out, so every relabelled cluster is auditable against its CellTypist
call and its scores.

The CellTypist input directory is a parameter (--celltypist-dir). Resumability belongs to the
step, through its .done markers.

Usage: cell_labels.py --dataset DS --celltypist-dir DIR --outdir DIR [--sample S]
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flexio
from celltype_annotate import ORGANISM, build_gene_adata, obs_table


PROGRAMS = {
    "erythroid": {
        "mouse": ["Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Alas2", "Gypa", "Slc4a1", "Car2"],
        "human": ["HBA1", "HBA2", "HBB", "ALAS2", "GYPA", "SLC4A1", "CA2"]},
    "plasmablast": {
        "mouse": ["Jchain", "Xbp1", "Prdm1", "Sdc1", "Derl3", "Mzb1", "Tnfrsf17"],
        "human": ["JCHAIN", "XBP1", "PRDM1", "SDC1", "DERL3", "MZB1", "TNFRSF17"]},
    "germinal_centre_B": {
        "mouse": ["Aicda", "Gcsam", "Fas", "Mef2b", "Rgs13", "S1pr2"],
        "human": ["AICDA", "GCSAM", "FAS", "MEF2B", "RGS13", "S1PR2"]},
    "proliferating": {
        "mouse": ["Mki67", "Top2a", "Birc5", "Stmn1", "Ccnb1", "Ube2c"],
        "human": ["MKI67", "TOP2A", "BIRC5", "STMN1", "CCNB1", "UBE2C"]},
}

# A cluster is relabelled only when its program score is this many robust SDs (MAD-based)
# above the median cluster score. The threshold is set from the observed distribution, which
# is sharply bimodal: across the mouse samples the per-cluster z for a program that is absent
# tops out around 6 (3.3, 4.2, 4.6, 5.2, 6.3), while every cluster carrying the program scores
# 15-42, with one borderline case at 10.1. Anything from 8 to 15 gives the same calls; 10 sits
# in the gap. A threshold near 3 relabelled a third of all clusters, which is why it is not 3.
Z_THRESHOLD = 10.0

# Spleen is full of red blood cells and red-pulp macrophages phagocytose them, so ambient and
# engulfed haemoglobin transcripts put an erythroid signature on non-erythroid clusters. A
# genuine erythroid cluster also has a far higher fraction of its counts in haemoglobin genes,
# which the campaign already computed as pct_counts_hb, so that is required as well.
HB_FOLD_OVER_MEDIAN = 5.0
CLUSTER_KEY = "leiden_res1_0"


def score_programs(dataset, sample, obs):
    """-> (per-cluster program score frame, cluster->celltypist label)."""
    import scanpy as sc

    organism = ORGANISM.get(dataset, "human")
    A = build_gene_adata(dataset, sample, obs, organism)   # log1p-normalised, human symbols
    # build_gene_adata maps mouse->human symbols, so score with the human gene lists
    present = set(A.var_names)
    scores = {}
    for prog, sets in PROGRAMS.items():
        genes = [g for g in sets["human"] if g in present]
        if not genes:
            continue
        sc.tl.score_genes(A, genes, score_name=f"_score_{prog}", ctrl_size=50)
        scores[prog] = A.obs[f"_score_{prog}"].values
    if CLUSTER_KEY not in A.obs.columns:
        return None, None, {}
    clu = A.obs[CLUSTER_KEY].astype(str).values
    df = pd.DataFrame(scores, index=A.obs_names)
    df["cluster"] = clu
    per_cluster = df.groupby("cluster").mean()
    used = {p: [g for g in s["human"] if g in present] for p, s in PROGRAMS.items()}
    return per_cluster, pd.Series(clu, index=A.obs_names), used


def refine(dataset, sample, outdir, ct_dir):
    ct_path = os.path.join(ct_dir, f"{sample}.celltypist.tsv.gz")
    if not os.path.exists(ct_path):
        raise SystemExit(f"[labels] {dataset}/{sample}: no CellTypist output {ct_path}")
    ct = pd.read_csv(ct_path, sep="\t", index_col=0)
    obs = obs_table(dataset, sample)
    per_cluster, clu, used = score_programs(dataset, sample, obs)
    if per_cluster is None:
        raise SystemExit(f"[labels] {dataset}/{sample}: no {CLUSTER_KEY} in the obs table")

    # robust z of each cluster's score within each program
    z = per_cluster.copy()
    for p in per_cluster.columns:
        v = per_cluster[p].values
        mad = np.median(np.abs(v - np.median(v))) * 1.4826
        z[p] = (v - np.median(v)) / (mad if mad > 0 else np.inf)

    base = "majority_voting" if "majority_voting" in ct.columns else "predicted_labels"
    ct["celltypist"] = ct[base].astype(str)
    ct["cluster"] = clu.reindex(ct.index).astype(str)

    # A cluster can genuinely carry two programs at once -- proliferating germinal-centre B
    # cells are the clearest case here -- so every program over the threshold is kept rather
    # than arbitrarily taking the argmax.
    hb = None
    if obs is not None and "pct_counts_hb" in obs.columns:
        o = obs.copy()
        o["_c"] = o[CLUSTER_KEY].astype(str)
        hb = o.groupby("_c")["pct_counts_hb"].mean()
    hb_median = float(np.median(hb.values)) if hb is not None and len(hb) else None

    remap = {}
    for c in z.index:
        called = [p for p in z.columns if z.loc[c, p] >= Z_THRESHOLD]
        if "erythroid" in called and hb is not None and hb_median is not None:
            enough_hb = hb.get(c, 0.0) >= max(hb_median * HB_FOLD_OVER_MEDIAN, 1e-9)
            if not enough_hb:
                called.remove("erythroid")
        remap[c] = "+".join(sorted(called)) if called else None
    ct["marker_program"] = ct["cluster"].map(lambda c: remap.get(c) or "")
    ct["label"] = np.where(ct["marker_program"] != "", ct["marker_program"], ct["celltypist"])

    os.makedirs(outdir, exist_ok=True)
    ct.to_csv(os.path.join(outdir, f"{sample}.labels.tsv.gz"), sep="\t")
    z.to_csv(os.path.join(outdir, f"{sample}.cluster_program_z.tsv"), sep="\t")
    relab = {c: v for c, v in remap.items() if v}
    print(f"[labels] {dataset}/{sample}: {len(ct)} cells, "
          f"{ct.label.nunique()} labels; {len(relab)} of {len(z)} clusters relabelled "
          f"{relab if relab else ''}", flush=True)
    return {"n_cells": int(len(ct)), "n_labels": int(ct.label.nunique()),
            "relabelled_clusters": relab, "markers_used": used,
            "label_counts": ct.label.value_counts().to_dict()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--sample")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--celltypist-dir", required=True)
    a = ap.parse_args()
    samples = [a.sample] if a.sample else [s for s, _ in flexio.samples_of(a.dataset)]
    res = {}
    for s in samples:
        res[s] = refine(a.dataset, s, a.outdir, a.celltypist_dir)
    with open(os.path.join(a.outdir, f"labels{'.' + a.sample if a.sample else ''}.json"), "w") as f:
        json.dump(res, f, indent=2)


if __name__ == "__main__":
    main()
