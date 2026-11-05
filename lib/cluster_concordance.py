#!/usr/bin/env python3
"""
Agreement between the clusterings that the sc-best-practices pipeline produced from each
tool's counts for the same Flex sample.

Five partition-similarity measures (ARI and NMI plus three more), computed on the cells
that survived QC in BOTH tools, for every Leiden resolution the pipeline ran:

  ARI  adjusted Rand index          — pair counting, adjusted for chance      (1 = identical)
  NMI  normalized mutual information — information shared, arithmetic mean normalizer
  AMI  adjusted mutual information   — NMI adjusted for chance (fairer at unequal k)
  V    V-measure                     — harmonic mean of homogeneity & completeness
       (homogeneity and completeness are reported separately as well)
  FMI  Fowlkes-Mallows index         — geometric mean of pair precision & recall
  NVI  normalized variation of information — a true metric; 0 = identical, 1 = independent

ARI, NMI, AMI, V and NVI are symmetric. Homogeneity and completeness are not: in every
pair the FIRST tool (alphabetical, so `cellranger` whenever it is present) is
`labels_true`, i.e. completeness < 1 means that tool's clusters got split by the other.

The pair also reports how the two cell sets themselves differ (Jaccard, per-tool QC
survivors, cluster counts), because a clustering can only agree on cells both tools kept.

Usage:
  cluster_concordance.py --indir DIR --dataset DS --sample S [--gene-axis common] \
      --out cluster_cmp.json
"""
import argparse
import glob
import itertools
import json
import os

import numpy as np
import pandas as pd
from sklearn.metrics import (adjusted_mutual_info_score, adjusted_rand_score,
                             fowlkes_mallows_score,
                             homogeneity_completeness_v_measure, mutual_info_score,
                             normalized_mutual_info_score)


def _entropy(labels):
    _, counts = np.unique(labels, return_counts=True)
    p = counts / counts.sum()
    return float(-(p * np.log(p)).sum())


def partition_metrics(a, b):
    """All five measures for two labelings of the same cells (plus components)."""
    a = np.asarray(a)
    b = np.asarray(b)
    hom, com, vms = homogeneity_completeness_v_measure(a, b)
    mi = mutual_info_score(a, b)
    ha, hb = _entropy(a), _entropy(b)
    vi = ha + hb - 2 * mi                       # variation of information (nats)
    joint = ha + hb - mi                        # H(A,B)
    return {
        "ARI": float(adjusted_rand_score(a, b)),
        "NMI": float(normalized_mutual_info_score(a, b, average_method="arithmetic")),
        "AMI": float(adjusted_mutual_info_score(a, b)),
        "V_measure": float(vms),
        "FMI": float(fowlkes_mallows_score(a, b)),
        "NVI": float(vi / joint) if joint > 0 else 0.0,
        "homogeneity": float(hom), "completeness": float(com),
        "VI_nats": float(vi),
        "n_clusters_a": int(len(np.unique(a))), "n_clusters_b": int(len(np.unique(b))),
        "n_cells": int(len(a)),
    }


def load_obs(indir, dataset, sample, gene_axis):
    out = {}
    pat = os.path.join(indir, f"{dataset}__{sample}__*__{gene_axis}.obs.tsv.gz")
    for p in sorted(glob.glob(pat)):
        tool = os.path.basename(p).split("__")[2]
        out[tool] = pd.read_csv(p, sep="\t", index_col=0)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--gene-axis", default="common")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    obs = load_obs(args.indir, args.dataset, args.sample, args.gene_axis)
    if len(obs) < 2:
        raise SystemExit(f"need >=2 tools, found {sorted(obs)} in {args.indir}")

    res_keys = sorted({c for df in obs.values() for c in df.columns
                       if c.startswith("leiden_res")})
    result = {
        "dataset": args.dataset, "sample": args.sample, "gene_axis": args.gene_axis,
        "tools": {t: {"n_cells_after_qc": int(len(df)),
                      "n_clusters": {k: int(df[k].nunique()) for k in res_keys
                                     if k in df},
                      "pct_doublets": (float(100 * (df["scDblFinder_class"] == "doublet").mean())
                                       if "scDblFinder_class" in df else None)}
                  for t, df in obs.items()},
        "pairs": {},
    }
    for t1, t2 in itertools.combinations(sorted(obs), 2):
        d1, d2 = obs[t1], obs[t2]
        common = d1.index.intersection(d2.index)
        pair = {
            "n_common_cells": int(len(common)),
            "cell_jaccard": float(len(common) / len(d1.index.union(d2.index))),
            "resolutions": {},
        }
        for k in res_keys:
            if k in d1 and k in d2 and len(common):
                pair["resolutions"][k] = partition_metrics(
                    d1.loc[common, k].to_numpy(), d2.loc[common, k].to_numpy())
        result["pairs"][f"{t1}_vs_{t2}"] = pair

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    json.dump(result, open(args.out, "w"), indent=2)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
