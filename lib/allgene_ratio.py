#!/usr/bin/env python3
"""Unspliced-to-spliced ratio per sample x cell label from the 10x `region` column.

Feeds Figure 2 panel c. For each sample and each
cell label, the probe counts of every gene that carries probes of BOTH kinds (`region_csv`
spliced and unspliced, as the vendor CSV labels them) are summed by kind. This is the vendor
annotation, not the 50-mer relocation (`genomic_class`) that the pre-mRNA instrument uses.

Outputs
  --out TSV    donor, cluster, spliced, unspliced, cells (one row per sample x label)
  --json JSON  pooled ratio per label (labels with >= --min-label-umi on the paired probes),
               per-sample ratios of the --order labels and whether their order holds per sample

Usage:
  allgene_ratio.py --probe-root DIR --label-root DIR --dataset DS --out TSV --json JSON
                   --order "Naive B cells" "germinal_centre_B+proliferating" "proliferating"
                   [--samples S ...] [--min-label-umi 5e5] [--min-point-umi 2e4]
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from probe_instruments import probe_paths, read_labels, write_doublet_log


def ratio_table(probe_root, label_root, dataset, samples=None, drop_doublets=False):
    import anndata as ad
    out = []
    for sample, path in probe_paths(probe_root, dataset, samples):
        a = ad.read_h5ad(path)
        lab = read_labels(label_root, dataset, sample, drop_doublets)
        common = [b for b in a.obs_names if b in lab.index]
        if not common:
            raise SystemExit(f"{sample}: no barcode shared between matrix and labels")
        a = a[common]
        L = lab.loc[common].values
        reg = a.var["region_csv"].astype(str).values
        bad = set(reg) - {"spliced", "unspliced"}
        if bad:
            raise SystemExit(f"{sample}: unexpected region_csv values {bad}")
        gid = a.var["gene_id"].astype(str).values
        kinds = pd.DataFrame({"g": gid, "r": reg}).groupby("g").r.nunique()
        both = set(kinds[kinds == 2].index)
        keep = np.array([g in both for g in gid])
        X = a.X.tocsr()
        sp_cols, un_cols = keep & (reg == "spliced"), keep & (reg == "unspliced")
        for c in pd.unique(L):
            m = L == c
            s = X[m]
            out.append((sample, c, float(s[:, sp_cols].sum()), float(s[:, un_cols].sum()),
                        int(m.sum())))
        print(f"[ratio] {sample}: {len(common)} labelled cells, {len(both)} genes with both "
              f"kinds", flush=True)
    return pd.DataFrame(out, columns=["donor", "cluster", "spliced", "unspliced", "cells"])


def summarise(t, order, min_label_umi, min_point_umi):
    P = t.groupby("cluster")[["spliced", "unspliced"]].sum()
    P["total"] = P.spliced + P.unspliced
    P["pooled_ratio"] = P.unspliced / P.spliced
    shown = P[P.total >= min_label_umi].sort_values("pooled_ratio", ascending=False)
    miss = [c for c in order if c not in P.index]
    if miss:
        raise SystemExit(f"order labels absent from the data: {miss}")
    per = {}
    for d, g in t.groupby("donor"):
        g = g.set_index("cluster")
        ok = (g.spliced > 0) & (g.spliced + g.unspliced >= min_point_umi)
        r = (g.unspliced / g.spliced)[ok]
        per[d] = {c: (float(r[c]) if c in r.index else None) for c in order}
    holds = {d: all(v[a] is not None and v[b] is not None and v[a] > v[b]
                    for a, b in zip(order, order[1:])) for d, v in per.items()}
    return {
        "pooled_ratio": {c: float(P.loc[c, "pooled_ratio"]) for c in shown.index},
        "pooled_total_umi": {c: float(P.loc[c, "total"]) for c in shown.index},
        "order": order,
        "order_pooled": {c: float(P.loc[c, "pooled_ratio"]) for c in order},
        "per_sample": per,
        "order_holds_per_sample": holds,
        "n_samples_order_holds": int(sum(holds.values())),
        "n_samples": len(holds),
        "min_label_umi": min_label_umi, "min_point_umi": min_point_umi,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe-root", required=True); ap.add_argument("--label-root", required=True)
    ap.add_argument("--dataset", required=True); ap.add_argument("--samples", nargs="+")
    ap.add_argument("--drop-doublets", action="store_true")
    ap.add_argument("--out", required=True); ap.add_argument("--json", required=True)
    ap.add_argument("--order", nargs="+", required=True)
    ap.add_argument("--min-label-umi", type=float, default=5e5)
    ap.add_argument("--min-point-umi", type=float, default=2e4)
    a = ap.parse_args()
    t = ratio_table(a.probe_root, a.label_root, a.dataset, a.samples, a.drop_doublets)
    if a.drop_doublets:
        write_doublet_log(a.out.replace(".tsv", ".doublets_removed.tsv"))
    t.to_csv(a.out, sep="\t", index=False)
    s = summarise(t, a.order, a.min_label_umi, a.min_point_umi)
    s["dataset"] = a.dataset
    json.dump(s, open(a.json, "w"), indent=1)
    for c, r in s["pooled_ratio"].items():
        print(f"  {c:40s} pooled {r:.4f}")
    print(f"[ratio] order {' > '.join(a.order)} holds in {s['n_samples_order_holds']}/"
          f"{s['n_samples']} samples")


if __name__ == "__main__":
    main()
