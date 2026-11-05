#!/usr/bin/env python3
"""Step 09 tables: Supp. Table 4, the 16plex Results numbers, sweep, seed control, equal-k exceptions.

Reads the per-sample JSONs under ROOT/{kmatch,kmatch_k*,deg,objects} and writes into OUT:
  matched_k_status.tsv            per sample: target k, each tool's k / resolution / search method
  exceptions.tsv                  every (sample, tool, k-run) whose k != target (header only when none)
  clustering_matched_k.tsv        per sample x pair: ARI NMI AMI V FMI NVI homogeneity completeness
  clustering_matched_k_sweep.tsv  the same at each requested k (+ _summary: means per k x pair)
  clustering_reseed_control.tsv   seed-control vs seed-0 Cell Ranger at the matched target k
  seedctl_by_resolution.tsv       the same at the fixed resolutions 0.25/0.5/1.0
  deg_overlap_matched_k.tsv       per sample x pair marker Jaccard / shared-of-top-N
  tbl_biol.tsv                    Supp. Table 4: per dataset, means over ALL 3 pairs (the published table's basis),
                                  plus the same over the 2 Cell Ranger pairs only (clustering and marker agreement)
  results_numbers.json            clustering and marker-agreement numbers for each dataset (16plex is the paper's)
Gates (exit 1): every expected sample has kmatch, deg, every sweep k and three objects; every expected
seed-control sample has a scored control.

Usage: assemble_biol.py --root DIR --expected TSV --expected-seedctl TSV --sweep-k "4 8 12" --top 50 --out DIR
"""
import argparse, csv, json, os, sys

import numpy as np

LABEL = {"Flexv2_16k_Human": "16k human PBMC v2 4plex", "16k_mouse_Spleno": "16k mouse splenocyte v2 4plex",
         "40k_DTC": "40k human DTC v1 4plex", "16plex_human": "128k human PBMC v1 16plex",
         "320k_GEMX": "320k human PBMC v1 16plex"}
ORDER = list(LABEL)
TOOLS = ("cellranger", "simpleaf", "cyto")
PAIRS = ("cellranger_vs_cyto", "cellranger_vs_simpleaf", "cyto_vs_simpleaf")
CR_PAIRS = PAIRS[:2]
MET = ["ARI", "NMI", "AMI", "V_measure", "FMI", "NVI", "homogeneity", "completeness", "VI_nats"]


def load(p):
    with open(p) as f:
        return json.load(f)


def write(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def mean(rows, col, pairs=None):
    v = [float(r[col]) for r in rows if (pairs is None or r["pair"] in pairs) and r.get(col) is not None]
    return float(np.mean(v)) if v else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--expected", required=True)
    ap.add_argument("--expected-seedctl", required=True)
    ap.add_argument("--sweep-k", required=True)
    ap.add_argument("--top", type=int, default=50)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    R = a.root
    exp = [tuple(l.rstrip("\n").split("\t")) for l in open(a.expected) if l.strip()]
    exp_sc = {tuple(l.rstrip("\n").split("\t")) for l in open(a.expected_seedctl) if l.strip()}
    sweep = a.sweep_k.split()
    shared = f"mean_shared_of_top{a.top}"

    bad = []
    for ds, s in exp:
        need = [f"kmatch/{ds}__{s}.json", f"deg/{ds}__{s}.json"] + \
               [f"kmatch_k{k}/{ds}__{s}.json" for k in sweep] + \
               [f"objects/{ds}__{s}__{t}__common.summary.json" for t in TOOLS]
        bad += [f"{ds}/{s}: missing {n}" for n in need if not os.path.exists(os.path.join(R, n))]
    for ds, s in sorted(exp_sc):
        p = os.path.join(R, "kmatch", f"{ds}__{s}.json")
        if os.path.exists(p) and not load(p).get("control"):
            bad.append(f"{ds}/{s}: seed control not scored in kmatch")
        if not os.path.exists(os.path.join(R, "objects", f"seedctl__{ds}__{s}__common.json")):
            bad.append(f"{ds}/{s}: no seedctl JSON")
    if bad:
        print("GATE FAIL:\n  " + "\n  ".join(bad))
        sys.exit(1)

    status, exc, km, deg, ctl, sw = [], [], [], [], [], []
    for ds, s in exp:
        K = load(os.path.join(R, "kmatch", f"{ds}__{s}.json"))
        D = load(os.path.join(R, "deg", f"{ds}__{s}.json"))
        summ = {t: load(os.path.join(R, "objects", f"{ds}__{s}__{t}__common.summary.json")) for t in TOOLS}
        row = {"dataset": ds, "sample": s, "target_k": K["target_k"]}
        for t in TOOLS:
            m = K["matched"][t]
            row.update({f"{t}_k": m["k"], f"{t}_resolution": m["resolution"], f"{t}_method": m["method"],
                        f"{t}_cells_after_qc": summ[t]["n_cells_after_qc"],
                        f"{t}_deg_clusters": D["n_clusters"].get(t)})
        row["equal_k"] = all(K["matched"][t]["k"] == K["target_k"] for t in TOOLS)
        status.append(row)
        for kdir, J in [("kmatch", K)] + [(f"kmatch_k{k}", load(os.path.join(R, f"kmatch_k{k}", f"{ds}__{s}.json")))
                                         for k in sweep]:
            for t, m in J["matched"].items():
                if m["k"] != J["target_k"]:
                    exc.append({"dataset": ds, "sample": s, "run": kdir, "tool": t, "target_k": J["target_k"],
                                "achieved_k": m["k"], "resolution": m["resolution"], "method": m["method"]})
            c = J.get("control")
            if kdir == "kmatch" and c and c["k"] != J["target_k"]:
                exc.append({"dataset": ds, "sample": s, "run": kdir, "tool": "cellranger_seedctl",
                            "target_k": J["target_k"], "achieved_k": c["k"], "resolution": c["resolution"],
                            "method": c["method"]})
            if kdir != "kmatch":
                for p, v in J["pairs"].items():
                    sw.append({"requested_k": J["target_k"], "dataset": ds, "sample": s, "pair": p,
                               "k_a": v["n_clusters_a"], "k_b": v["n_clusters_b"], "n_cells": v["n_cells"],
                               **{m_: v[m_] for m_ in MET}})
        for p, v in K["pairs"].items():
            km.append({"dataset": ds, "sample": s, "target_k": K["target_k"], "pair": p,
                       "k_a": v["n_clusters_a"], "k_b": v["n_clusters_b"], "n_cells": v["n_cells"],
                       **{m_: v[m_] for m_ in MET}})
        if K.get("control"):
            c = K["control"]
            ctl.append({"dataset": ds, "sample": s, "target_k": K["target_k"], "k": c["k"],
                        "resolution": c["resolution"], "n_cells": c["n_cells"], **{m_: c[m_] for m_ in MET}})
        for p, v in D["pairs"].items():
            deg.append({"dataset": ds, "sample": s, "target_k": D["target_k"], "pair": p,
                        "matched_clusters": v["matched_clusters"],
                        "mean_marker_jaccard": v["mean_marker_jaccard"],
                        "median_marker_jaccard": v["median_marker_jaccard"], shared: v[shared]})

    sc_res = []
    for ds, s in sorted(exp_sc):
        J = load(os.path.join(R, "objects", f"seedctl__{ds}__{s}__common.json"))
        for res, v in J["resolutions"].items():
            sc_res.append({"dataset": ds, "sample": s, "resolution": res, "n_cells": J["n_common_cells"],
                           **{m_: v[m_] for m_ in MET}})

    scols = ["dataset", "sample", "target_k", "equal_k"] + \
            [f"{t}_{c}" for t in TOOLS for c in ("k", "resolution", "method", "cells_after_qc", "deg_clusters")]
    write(os.path.join(a.out, "matched_k_status.tsv"), status, scols)
    write(os.path.join(a.out, "exceptions.tsv"), exc,
          ["dataset", "sample", "run", "tool", "target_k", "achieved_k", "resolution", "method"])
    pc = ["dataset", "sample", "target_k", "pair", "k_a", "k_b", "n_cells"] + MET
    write(os.path.join(a.out, "clustering_matched_k.tsv"), km, pc)
    write(os.path.join(a.out, "clustering_matched_k_sweep.tsv"), sw, ["requested_k"] + pc[:2] + pc[3:])
    summ_rows = []
    for k in sweep:
        for p in PAIRS:
            rs = [r for r in sw if str(r["requested_k"]) == k and r["pair"] == p]
            summ_rows.append({"requested_k": k, "pair": p, "n_samples": len(rs),
                              "n_equal_k": sum(r["k_a"] == r["k_b"] == int(k) for r in rs),
                              **{m_: mean(rs, m_) for m_ in MET}})
    write(os.path.join(a.out, "clustering_matched_k_sweep_summary.tsv"), summ_rows,
          ["requested_k", "pair", "n_samples", "n_equal_k"] + MET)
    write(os.path.join(a.out, "clustering_reseed_control.tsv"), ctl,
          ["dataset", "sample", "target_k", "k", "resolution", "n_cells"] + MET)
    write(os.path.join(a.out, "seedctl_by_resolution.tsv"), sc_res,
          ["dataset", "sample", "resolution", "n_cells"] + MET)
    write(os.path.join(a.out, "deg_overlap_matched_k.tsv"), deg,
          ["dataset", "sample", "target_k", "pair", "matched_clusters", "mean_marker_jaccard",
           "median_marker_jaccard", shared])

    tbl, nums = [], {}
    for ds in [d for d in ORDER if any(x[0] == d for x in exp)]:
        k_ = [r for r in km if r["dataset"] == ds]
        d_ = [r for r in deg if r["dataset"] == ds]
        st = [r for r in status if r["dataset"] == ds]
        tbl.append({"dataset": ds, "label": LABEL[ds], "n_samples": len(st),
                    "n_pairs": len(k_), "n_equal_k_samples": sum(r["equal_k"] for r in st),
                    "homogeneity": mean(k_, "homogeneity"), "completeness": mean(k_, "completeness"),
                    "ARI": mean(k_, "ARI"),
                    "marker_jaccard": mean(d_, "mean_marker_jaccard"), "shared_top": mean(d_, shared),
                    "ARI_cr_pairs": mean(k_, "ARI", CR_PAIRS),
                    "homogeneity_cr_pairs": mean(k_, "homogeneity", CR_PAIRS),
                    "completeness_cr_pairs": mean(k_, "completeness", CR_PAIRS),
                    "marker_jaccard_cr_pairs": mean(d_, "mean_marker_jaccard", CR_PAIRS),
                    "shared_top_cr_pairs": mean(d_, shared, CR_PAIRS)})
        cr_cells = [r["cellranger_cells_after_qc"] for r in st]
        nums[ds] = {
            "n_samples": len(st), "cr_cells_after_qc_min": min(cr_cells), "cr_cells_after_qc_max": max(cr_cells),
            "mean_target_k": float(np.mean([r["target_k"] for r in st])),
            "n_pairwise_comparisons": len(k_),
            # equal k = both tools' clusterings reached the target (matched k on each tool's own cells),
            # not k_a == k_b, which counts clusters on the shared cells only
            "n_pairs_equal_k": sum(all(m[f"{t}_k"] == m["target_k"] for t in r["pair"].split("_vs_"))
                                   for r in k_ for m in st if m["sample"] == r["sample"]),
            "homogeneity_all_pairs": mean(k_, "homogeneity"), "completeness_all_pairs": mean(k_, "completeness"),
            "shared_top_cr_pairs": mean(d_, shared, CR_PAIRS), "shared_top_all_pairs": mean(d_, shared),
            "marker_jaccard_cr_pairs": mean(d_, "mean_marker_jaccard", CR_PAIRS),
            "marker_jaccard_all_pairs": mean(d_, "mean_marker_jaccard"),
            "marker_jaccard_by_pair": {p: mean(d_, "mean_marker_jaccard", (p,)) for p in PAIRS},
            "seed_control_samples": sum(r["dataset"] == ds for r in ctl),
            "seed_control_ARI_mean": mean([r for r in ctl if r["dataset"] == ds], "ARI")}
    tcols = ["dataset", "label", "n_samples", "n_pairs", "n_equal_k_samples", "homogeneity", "completeness",
             "ARI", "marker_jaccard", "shared_top", "homogeneity_cr_pairs", "completeness_cr_pairs",
             "ARI_cr_pairs", "marker_jaccard_cr_pairs", "shared_top_cr_pairs"]
    write(os.path.join(a.out, "tbl_biol.tsv"), tbl, tcols)
    nums["_definitions"] = {
        "tbl_biol": "per-dataset unweighted means over every (sample, pair) row; *_cr_pairs = Cell Ranger pairs only",
        "homogeneity/completeness": "first tool of the pair name (alphabetical) is labels_true",
        "seed_control_samples_total": len(ctl), "equal_k_exceptions": len(exc)}
    json.dump(nums, open(os.path.join(a.out, "results_numbers.json"), "w"), indent=2)

    for r in tbl:
        print(f"{r['label']:32s} hom={r['homogeneity']:.3f} com={r['completeness']:.3f} ARI={r['ARI']:.3f} "
              f"jac={r['marker_jaccard']:.3f} shared={r['shared_top']:.1f} | CR-pairs jac="
              f"{r['marker_jaccard_cr_pairs']:.3f} shared={r['shared_top_cr_pairs']:.1f} "
              f"| equal-k {r['n_equal_k_samples']}/{r['n_samples']}")
    print(f"equal-k exceptions: {len(exc)}" + "".join(
        f"\n  EXCEPTION {e['run']} {e['dataset']}/{e['sample']} {e['tool']}: k={e['achieved_k']} target "
        f"{e['target_k']} ({e['method']})" for e in exc))


if __name__ == "__main__":
    main()
