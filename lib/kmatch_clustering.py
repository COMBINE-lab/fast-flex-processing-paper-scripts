#!/usr/bin/env python3
"""
Cluster-count-matched agreement: compare partitions at the SAME number of clusters.

match_k does not give up when bisection misses the target. Leiden's cluster count is a
NON-monotone step function
of resolution, so bisection can step over every window where k == target and settle on target +- 1
(16plex PBMC_15 cyto: bisection parked at r=0.931/k=21 while k=20 holds at r=0.825, 0.85, 0.885, ...).
The search is now:
  1. bisection over [lo, hi] exactly as before (so every sample it already solved is unchanged);
  2. on a miss, a dense grid (scan_step over [lo, hi], plus a log-spaced grid from scan_min up to lo
     for coarse targets such as k=4) visited outward from bisection's last midpoint, first hit wins;
  3. still no hit: a fine scan (fine_step) inside the grid intervals whose ends straddle the target,
     nearest first;
  4. still no hit: the closest k is kept and the tool is flagged `exception` in the JSON -- never silent.
The Leiden seed stays fixed (random_state=0, the paper's seed): deg_overlap.py re-runs Leiden at the
stored resolution with that seed, so a match found under another seed would not reproduce there.

Each tool's saved .h5ad already holds the PCA and neighbour graph, so nothing upstream is recomputed.
Target k defaults to Cell Ranger's k at resolution 1.0. The seed-1 control object, where it exists, is
matched to the same k.

Usage: kmatch_clustering.py --indir DIR --dataset DS --sample S --out JSON [--target-k K]
           [--gene-axis common] [--seed 1] [--lo .02 --hi 6 --bisect-steps 18
            --scan-step .005 --scan-min 1e-4 --fine-step 1e-4]
"""
import argparse
import itertools
import json
import os
import sys

import numpy as np
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cluster_concordance import partition_metrics

SEED = 0
TOOLS = ("cellranger", "simpleaf", "cyto")


def leiden_k(adata, res, key):
    sc.tl.leiden(adata, key_added=key, resolution=res, flavor="igraph",
                 n_iterations=2, random_state=SEED)
    return adata.obs[key].nunique()


def match_k(adata, target, lo=0.02, hi=6.0, tries=18, scan_step=0.005, scan_min=1e-4,
            fine_step=1e-4, max_fine_intervals=20):
    """-> (labels, resolution, k, method, n_leiden_runs); method in bisection|dense_scan|fine_scan|unreached."""
    seen = {}

    def run(r):
        r = float(r)
        if r not in seen:
            seen[r] = (leiden_k(adata, r, "_km"), adata.obs["_km"].to_numpy())
        return seen[r][0]

    best, mid, lo0, hi0 = None, None, lo, hi
    for _ in range(tries):                       # 1. bisection, identical to the benchmark version
        mid = (lo + hi) / 2
        k = run(mid)
        cand = (abs(k - target), abs(mid - 1.0), mid)
        if best is None or cand[:2] < best[:2]:
            best = cand
        if k == target:
            return seen[mid][1], mid, k, "bisection", len(seen)
        if k < target:
            lo = mid
        else:
            hi = mid

    def hit(points, method):
        for r in points:
            if run(r) == target:
                return seen[float(r)][1], float(r), target, method, len(seen)
        return None

    center = mid                                  # 2. dense grid, outward from the bracket
    grid = np.round(np.arange(lo0, hi0 + scan_step / 2, scan_step), 6)
    low = np.round(np.geomspace(scan_min, lo0, 40, endpoint=False), 6)
    pts = sorted(set(grid.tolist()) | set(low.tolist()), key=lambda r: (abs(r - center), r))
    got = hit(pts, "dense_scan")
    if got:
        return got

    rs = sorted(seen)                             # 3. fine scan where the grid straddles the target
    straddle = [(a, b) for a, b in zip(rs, rs[1:])
                if (seen[a][0] - target) * (seen[b][0] - target) < 0]
    straddle.sort(key=lambda ab: abs((ab[0] + ab[1]) / 2 - center))
    for a, b in straddle[:max_fine_intervals]:
        got = hit(np.round(np.arange(a + fine_step, b, fine_step), 6), "fine_scan")
        if got:
            return got

    r = min(seen, key=lambda r: (abs(seen[r][0] - target), abs(r - 1.0)))   # 4. closest, flagged
    return seen[r][1], r, seen[r][0], "unreached", len(seen)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--gene-axis", default="common")
    ap.add_argument("--seed", type=int, default=1, help="seed of the control object to score")
    ap.add_argument("--out", required=True)
    ap.add_argument("--target-k", type=int, default=None,
                    help="cluster count to match every tool to (default: Cell Ranger's k at res 1.0)")
    ap.add_argument("--lo", type=float, default=0.02)
    ap.add_argument("--hi", type=float, default=6.0)
    ap.add_argument("--bisect-steps", type=int, default=18)
    ap.add_argument("--scan-step", type=float, default=0.005)
    ap.add_argument("--scan-min", type=float, default=1e-4)
    ap.add_argument("--fine-step", type=float, default=1e-4)
    args = ap.parse_args()
    kw = dict(lo=args.lo, hi=args.hi, tries=args.bisect_steps, scan_step=args.scan_step,
              scan_min=args.scan_min, fine_step=args.fine_step)

    def path(tool, seed=None):
        t = f"{args.dataset}__{args.sample}__{tool}__{args.gene_axis}"
        if seed:
            t += f"__seed{seed}"
        return os.path.join(args.indir, t + ".h5ad")

    objs = {t: sc.read_h5ad(path(t)) for t in TOOLS if os.path.exists(path(t))}
    if set(objs) != set(TOOLS):
        raise SystemExit(f"need all of {TOOLS} for {args.dataset}/{args.sample}, have {sorted(objs)}")

    target = args.target_k or int(objs["cellranger"].obs["leiden_res1_0"].nunique())
    res = {"dataset": args.dataset, "sample": args.sample, "gene_axis": args.gene_axis,
           "target_k": target, "search": kw, "matched": {}, "pairs": {}, "control": None,
           "exceptions": []}

    labels = {}
    for tool, ad in objs.items():
        lab, r, k, method, n = match_k(ad, target, **kw)
        labels[tool] = (ad.obs_names.to_numpy(), lab)
        res["matched"][tool] = {"resolution": float(r), "k": int(k), "method": method,
                                "leiden_runs": n,
                                "k_at_res1": int(ad.obs["leiden_res1_0"].nunique())}
        if k != target:
            res["exceptions"].append(tool)

    for a, b in itertools.combinations(sorted(labels), 2):
        na, la = labels[a]
        nb, lb = labels[b]
        common = np.intersect1d(na, nb)
        ia = {n: i for i, n in enumerate(na)}
        ib = {n: i for i, n in enumerate(nb)}
        res["pairs"][f"{a}_vs_{b}"] = partition_metrics(
            la[[ia[c] for c in common]], lb[[ib[c] for c in common]])

    ctl = path("cellranger", args.seed)
    if os.path.exists(ctl):
        ad = sc.read_h5ad(ctl)
        lab, r, k, method, n = match_k(ad, target, **kw)
        na, la = labels["cellranger"]
        common = np.intersect1d(na, ad.obs_names.to_numpy())
        ia = {n_: i for i, n_ in enumerate(na)}
        ib = {n_: i for i, n_ in enumerate(ad.obs_names.to_numpy())}
        res["control"] = {"resolution": float(r), "k": int(k), "method": method,
                          **partition_metrics(la[[ia[c] for c in common]],
                                              lab[[ib[c] for c in common]])}
        if k != target:
            res["exceptions"].append(f"cellranger__seed{args.seed}")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    json.dump(res, open(args.out, "w"), indent=2)
    line = " ".join(f"{t}:k={v['k']}@r={v['resolution']:.4f}({v['method']})"
                    for t, v in res["matched"].items())
    print(f"{args.dataset}/{args.sample} target_k={target} {line}")
    for p, m in res["pairs"].items():
        print(f"   {p:26s} ARI={m['ARI']:.3f} NMI={m['NMI']:.3f} hom={m['homogeneity']:.3f} "
              f"com={m['completeness']:.3f}")
    if res["exceptions"]:
        print(f"   EXCEPTION: target k={target} not reached for {res['exceptions']}")


if __name__ == "__main__":
    main()
