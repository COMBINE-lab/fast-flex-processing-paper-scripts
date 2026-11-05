#!/usr/bin/env python3
"""
Do the three tools yield the same marker genes?

Runs differential expression per cluster on each tool's object using the CLUSTER-COUNT-
MATCHED labels (so every tool has the same k), matches clusters between tools by maximum
overlap, and reports the overlap of each matched cluster's top-N markers.

Cluster labels are arbitrary integers, so a raw label-to-label comparison is meaningless:
clusters are paired with the Hungarian assignment on the contingency table (maximising
shared cells), and only then are their marker sets compared.

DE follows the book's annotation chapter: Wilcoxon rank-sum on the log-normalised layer.

Usage: deg_overlap.py --indir DIR --kmatch-dir DIR --dataset DS --sample S
                      [--top 50] --out JSON
"""
import argparse
import itertools
import json
import os
import sys

import numpy as np
import scanpy as sc
from scipy.optimize import linear_sum_assignment

TOOLS = ("cellranger", "simpleaf", "cyto")


def markers(adata, labels, top, layer="log1p_norm"):
    adata = adata.copy()
    adata.obs["_k"] = [str(x) for x in labels]
    if layer in adata.layers:
        adata.X = adata.layers[layer]
    keep = adata.obs["_k"].value_counts()
    keep = set(keep[keep >= 3].index)            # Wilcoxon needs a few cells per group
    adata = adata[adata.obs["_k"].isin(keep)].copy()
    sc.tl.rank_genes_groups(adata, "_k", method="wilcoxon")
    out = {}
    for g in adata.obs["_k"].cat.categories if hasattr(adata.obs["_k"], "cat") else sorted(keep):
        try:
            names = sc.get.rank_genes_groups_df(adata, group=str(g)).head(top)["names"]
            out[str(g)] = list(names)
        except (KeyError, ValueError):
            continue
    return out


def match_clusters(la, lb, names_a, names_b):
    """Hungarian pairing of cluster ids by shared cells on the common barcodes."""
    common = np.intersect1d(names_a, names_b)
    ia = {n: i for i, n in enumerate(names_a)}
    ib = {n: i for i, n in enumerate(names_b)}
    a = np.array([la[ia[c]] for c in common])
    b = np.array([lb[ib[c]] for c in common])
    ua, ub = sorted(set(a)), sorted(set(b))
    M = np.zeros((len(ua), len(ub)))
    for i, x in enumerate(ua):
        for j, y in enumerate(ub):
            M[i, j] = np.sum((a == x) & (b == y))
    r, c = linear_sum_assignment(-M)
    return [(str(ua[i]), str(ub[j]), int(M[i, j])) for i, j in zip(r, c)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--kmatch-dir", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--gene-axis", default="common")
    ap.add_argument("--top", type=int, default=50)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    km = os.path.join(args.kmatch_dir, f"{args.dataset}__{args.sample}.json")
    if not os.path.exists(km):
        raise SystemExit(f"no k-matched labels for {args.dataset}/{args.sample}")
    kinfo = json.load(open(km))

    objs, labels, obsnames, mk = {}, {}, {}, {}
    for tool in TOOLS:
        p = os.path.join(args.indir,
                         f"{args.dataset}__{args.sample}__{tool}__{args.gene_axis}.h5ad")
        if not os.path.exists(p) or tool not in kinfo["matched"]:
            continue
        ad = sc.read_h5ad(p)
        res = kinfo["matched"][tool]["resolution"]
        key = "_km"
        sc.tl.leiden(ad, key_added=key, resolution=res, flavor="igraph",
                     n_iterations=2, random_state=0)
        objs[tool] = ad
        labels[tool] = ad.obs[key].to_numpy()
        obsnames[tool] = ad.obs_names.to_numpy()
        mk[tool] = markers(ad, labels[tool], args.top)

    res = {"dataset": args.dataset, "sample": args.sample, "target_k": kinfo["target_k"],
           "top": args.top, "n_clusters": {t: len(m) for t, m in mk.items()}, "pairs": {}}
    for a, b in itertools.combinations(sorted(mk), 2):
        pairs = match_clusters(labels[a], labels[b], obsnames[a], obsnames[b])
        jac, shared = [], []
        for ca, cb, n in pairs:
            A, B = set(mk[a].get(ca, [])), set(mk[b].get(cb, []))
            if not A or not B:
                continue
            jac.append(len(A & B) / len(A | B))
            shared.append(len(A & B))
        res["pairs"][f"{a}_vs_{b}"] = {
            "matched_clusters": len(jac),
            "mean_marker_jaccard": float(np.mean(jac)) if jac else None,
            "median_marker_jaccard": float(np.median(jac)) if jac else None,
            f"mean_shared_of_top{args.top}": float(np.mean(shared)) if shared else None,
        }
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    json.dump(res, open(args.out, "w"), indent=2)
    print(f"{args.dataset}/{args.sample} k={res['target_k']}")
    for p, m in res["pairs"].items():
        print(f"   {p:26s} clusters={m['matched_clusters']:3d} "
              f"marker Jaccard mean={m['mean_marker_jaccard']} "
              f"shared/{args.top}={m[f'mean_shared_of_top{args.top}']}")


if __name__ == "__main__":
    main()
