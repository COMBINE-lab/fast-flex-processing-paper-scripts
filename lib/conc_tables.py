#!/usr/bin/env python3
"""Step 08 tables from the per-sample compare_consistency.py JSONs.

Writes into OUTDIR:
  per_sample.tsv                      one row per (sample, pair): the full alevin-fry-paper metric set
  table_s1_concordance_detail.tsv     per-library unweighted means,
                                      then ALL = unweighted mean over every sample
  pair_definitions.tsv                which tool is reference and which is test in each pair (UMI ratio = test/ref)
  text_ranges.json                    the count-concordance ranges (per-sample min/max over all pair rows) and
                                      the L513 mouse gene-Spearman in both conventions
  cr_provenance.tsv                   the Cell Ranger matrix behind every sample: path, software_version
                                      read from the h5 attributes, sha256

Fails (exit 1) unless every expected (dataset, sample) has a JSON with all three pairs and no error.

Usage: conc_tables.py --jsons DIR --expected TSV(dataset<TAB>sample) --out OUTDIR
"""
import argparse, csv, glob, hashlib, json, os, sys

import h5py
import numpy as np

PAIRS = [("cellranger_vs_simpleaf", "Cell Ranger vs simpleaf", "cellranger", "simpleaf"),
         ("cellranger_vs_cyto", "Cell Ranger vs cyto", "cellranger", "cyto"),
         ("simpleaf_vs_cyto", "simpleaf vs cyto", "simpleaf", "cyto")]
LIBS = ["Flexv2_16k_Human", "16k_mouse_Spleno", "40k_DTC", "16plex_human", "320k_GEMX"]
METRICS = ["n_common_cells", "n_common_genes", "mean_cell_spearman", "median_cell_spearman",
           "MARD_dropNA", "MARD_NA0", "mean_rFP_per_cell", "mean_rFN_per_cell",
           "gene_pearson_log1p", "gene_spearman", "count_concordance_pct",
           "total_umi_ref", "total_umi_test", "total_umi_ratio_test_over_ref",
           "cellranger_cells_recovered_pct"]
# Supp Table S1 column -> per-sample metric, decimals (the published file's rounding)
S1 = [("mean_cell_spearman", "mean_cell_spearman", 4), ("MARD_dropNA", "MARD_dropNA", 4),
      ("MARD_NA0", "MARD_NA0", 4), ("rFP_per_cell", "mean_rFP_per_cell", 4),
      ("rFN_per_cell", "mean_rFN_per_cell", 4), ("gene_pearson_log1p", "gene_pearson_log1p", 4),
      ("count_concordance_pct", "count_concordance_pct", 2),
      ("UMI_ratio_test_over_ref", "total_umi_ratio_test_over_ref", 4)]


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def cr_h5(path):
    if os.path.isfile(path):
        return path
    hit = sorted(glob.glob(os.path.join(path, "*sample_filtered_feature_bc_matrix.h5")))
    if len(hit) != 1:
        raise SystemExit(f"cannot pin one filtered h5 under {path}: {hit}")
    return hit[0]


def write(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jsons", required=True)
    ap.add_argument("--expected", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    expected = [tuple(l.rstrip("\n").split("\t")) for l in open(a.expected) if l.strip()]
    data, bad = {}, []
    for ds, s in expected:
        f = os.path.join(a.jsons, f"{ds}__{s}.json")
        if not os.path.exists(f):
            bad.append(f"{ds}/{s}: no JSON")
            continue
        d = json.load(open(f))
        p = d.get("pairs", {})
        miss = [k for k, *_ in PAIRS if k not in p or "error" in p[k]]
        if miss:
            bad.append(f"{ds}/{s}: missing/errored pairs {miss}")
        data[(ds, s)] = d
    if bad:
        print("GATE FAIL:\n  " + "\n  ".join(bad))
        sys.exit(1)

    # per-sample rows
    rows = []
    for (ds, s), d in data.items():
        for k, _, ref, test in PAIRS:
            m = d["pairs"][k]
            rows.append([ds, s, k, ref, test] + [m.get(c, "") for c in METRICS])
    write(os.path.join(a.out, "per_sample.tsv"), rows,
          ["dataset", "sample", "pair", "ref", "test"] + METRICS)
    write(os.path.join(a.out, "pair_definitions.tsv"),
          [[lab, k, ref, test, f"{test}/{ref}"] for k, lab, ref, test in PAIRS],
          ["s1_pair_label", "json_pair", "reference", "test", "UMI_ratio_test_over_ref"])

    # Table S1: per-library unweighted means, then ALL
    s1 = []
    libs = [l for l in LIBS if any(ds == l for ds, _ in data)]
    for lib in libs + ["ALL"]:
        keys = [k for k in data if lib == "ALL" or k[0] == lib]
        for k, lab, *_ in PAIRS:
            vals = [data[x]["pairs"][k] for x in keys]
            s1.append([lib, lab, len(vals)] +
                      [f"{np.mean([float(v[src]) for v in vals]):.{dp}f}" for _, src, dp in S1])
    write(os.path.join(a.out, "table_s1_concordance_detail.tsv"), s1,
          ["library", "pair", "n_samples"] + [c for c, _, _ in S1])

    # count concordance: per-sample min/max over every pair row; the mouse diagnostic
    def rng(c):
        v = [float(d["pairs"][k][c]) for d in data.values() for k, *_ in PAIRS]
        return {"min": min(v), "max": max(v), "n_rows": len(v)}
    mouse = [d["pairs"]["cellranger_vs_simpleaf"] for (ds, _), d in data.items()
             if ds == "16k_mouse_Spleno"]
    txt = {"n_samples": len(data),
           "mean_cell_spearman": rng("mean_cell_spearman"),
           "count_concordance_pct": rng("count_concordance_pct"),
           "MARD_NA0": rng("MARD_NA0"),
           "definition": "per-sample value (mean over CR-called cells), min/max over all samples x 3 pairs"}
    if mouse:
        txt["mouse_cellranger_vs_simpleaf"] = {
            "mean_of_per_sample_count_concordance_pct": float(np.mean([m["count_concordance_pct"] for m in mouse])),
            "mean_of_per_sample_gene_spearman": float(np.mean([m["gene_spearman"] for m in mouse])),
            "n_samples": len(mouse),
            "note": "per-sample means; the pooled 16,372-cell convention is not computed here"}
    json.dump(txt, open(os.path.join(a.out, "text_ranges.json"), "w"), indent=2)

    # Cell Ranger provenance, read from the matrices themselves
    prov = []
    for (ds, s), d in data.items():
        h5 = cr_h5(d["cellranger_matrix"])
        with h5py.File(h5, "r") as h:
            ver = h.attrs.get("software_version", b"")
        ver = ver.decode() if isinstance(ver, bytes) else str(ver)
        prov.append([ds, s, os.path.realpath(h5), ver, sha256(h5)])
    write(os.path.join(a.out, "cr_provenance.tsv"), prov,
          ["dataset", "sample", "cellranger_filtered_h5", "software_version", "sha256"])
    vers = sorted({(r[0], r[3]) for r in prov})
    print("Cell Ranger software_version per dataset:", vers)
    print(f"tables -> {a.out}: {len(rows)} pair rows, {len(data)} samples")


if __name__ == "__main__":
    main()
