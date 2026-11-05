#!/usr/bin/env python3
"""
How much of a tool-vs-tool clustering disagreement is just the pipeline's own randomness?

Runs the same comparison as cluster_concordance.py, but between two runs of the SAME tool's
matrix that differ only in seed. Leiden, scDblFinder and the SoupX pre-clustering are all
stochastic, so this is the ceiling any tool pair can reach on that sample: an ARI of 0.78
between two tools means little if two seeds on one matrix also give 0.78.

Usage:
  seed_control.py --indir DIR --dataset DS --sample S [--tool cellranger]
                  [--gene-axis common] [--seed 1] --out control.json
"""
import argparse
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cluster_concordance import partition_metrics


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--tool", default="cellranger")
    ap.add_argument("--gene-axis", default="common")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    base = os.path.join(args.indir,
                        f"{args.dataset}__{args.sample}__{args.tool}__{args.gene_axis}")
    a_path = base + ".obs.tsv.gz"
    b_path = f"{base}__seed{args.seed}.obs.tsv.gz"
    for p in (a_path, b_path):
        if not os.path.exists(p):
            raise SystemExit(f"missing {p}")
    a = pd.read_csv(a_path, sep="\t", index_col=0)
    b = pd.read_csv(b_path, sep="\t", index_col=0)
    common = a.index.intersection(b.index)

    res = {"dataset": args.dataset, "sample": args.sample, "tool": args.tool,
           "gene_axis": args.gene_axis, "seeds": [0, args.seed],
           "n_common_cells": int(len(common)),
           "cell_jaccard": float(len(common) / len(a.index.union(b.index))),
           "resolutions": {}}
    for k in [c for c in a.columns if c.startswith("leiden_res")]:
        if k in b:
            res["resolutions"][k] = partition_metrics(a.loc[common, k].to_numpy(),
                                                      b.loc[common, k].to_numpy())
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    json.dump(res, open(args.out, "w"), indent=2)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
