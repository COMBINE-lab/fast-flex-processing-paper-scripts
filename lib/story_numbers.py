#!/usr/bin/env python3
"""The per-gene probe facts behind Figure 2 panels a and b (YBX3, ZBTB16, PHACTR1).

Reads the paired PyDESeq2 gene-vs-probe table of one label (pseudobulk_de.py) and the
pseudobulk counts, and writes, for every probe of every --genes gene: probe id, genomic class,
isoform set and whether it contains all --tx transcripts, UMI in the label / the rest / the whole
library, fold change, padj, whether it passes the minimum-UMI floor, and the call with and
without the floor; plus the gene fold change and padj, and the probe count of each gene.

Usage: story_numbers.py --masked TSV --pb DIR --label L --genes G ... --tx T ... --out JSON
"""
import argparse
import json
import os

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--masked", required=True); ap.add_argument("--pb", required=True)
    ap.add_argument("--label", required=True); ap.add_argument("--genes", nargs="+", required=True)
    ap.add_argument("--tx", nargs="+", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--min-probe-umi", type=int, required=True)
    o = ap.parse_args()

    d = pd.read_csv(o.masked, sep="\t", index_col=0)
    M = np.load(os.path.join(o.pb, "probe_counts.npy"))
    meta = pd.read_csv(os.path.join(o.pb, "units.tsv"), sep="\t")
    var = pd.read_csv(os.path.join(o.pb, "probes.tsv"), sep="\t", index_col=0)
    var["tx_set"] = var["tx_set"].fillna("")
    lib_umi = pd.Series(M.sum(axis=0), index=var.index)
    in_label = pd.Series(M[(meta.label == o.label).values].sum(axis=0), index=var.index)

    res = {"label": o.label, "min_probe_umi": o.min_probe_umi, "genes": {}}
    for g in o.genes:
        pv = var[var.gene_name == g]
        if pv.empty:
            raise SystemExit(f"{g}: not in the probe panel")
        rows = d[d.gene_name == g]
        probes = []
        for pid, v in pv.iterrows():
            r = rows.loc[pid] if pid in rows.index else None
            tx = [t for t in v.tx_set.split(",") if t]
            probes.append({
                "probe_id": pid, "probe": pid.split("|")[-1],
                "genomic_class": v.genomic_class, "region_csv": v.region_csv,
                "n_tx_hit": int(v.n_tx_hit), "n_tx_gene": int(v.n_tx_gene),
                "reads_all_tx": all(t in tx for t in o.tx),
                "umi_library": float(lib_umi[pid]), "umi_label_pseudobulk": float(in_label[pid]),
                "tested": r is not None,
                "umi_target": None if r is None else float(r.umi_target),
                "umi_rest": None if r is None else float(r.umi_rest),
                "fold_change": None if r is None or pd.isna(r.probe_log2FoldChange)
                else float(2 ** r.probe_log2FoldChange),
                "padj": None if r is None or pd.isna(r.probe_padj) else float(r.probe_padj),
                "passes_floor": None if r is None else bool(r.passes_floor),
                "sig_without_floor": None if r is None else bool(r.probe_sig_nofloor),
                "sig": None if r is None else bool(r.probe_sig),
                "masked": None if r is None else bool(r.masked)})
        probes.sort(key=lambda p: -p["umi_library"])
        g0 = rows.iloc[0] if len(rows) else None
        res["genes"][g] = {
            "n_probes": len(pv),
            "gene_fold_change": None if g0 is None else float(2 ** g0.gene_log2FoldChange),
            "gene_padj": None if g0 is None else float(g0.gene_padj),
            "n_probes_sig": sum(bool(p["sig"]) for p in probes),
            "n_probes_sig_without_floor": sum(bool(p["sig_without_floor"]) for p in probes),
            "probes": probes}
        print(f"{g}: {len(pv)} probes, gene x{res['genes'][g]['gene_fold_change']:.3f} "
              f"padj {res['genes'][g]['gene_padj']:.3g}")
        for p in probes:
            fc = "-" if p["fold_change"] is None else f"x{p['fold_change']:.3f}"
            pa = "-" if p["padj"] is None else f"{p['padj']:.3g}"
            print(f"   {p['probe']} {p['genomic_class']:11s} tx {p['n_tx_hit']}/{p['n_tx_gene']}"
                  f" reads_all_tx={p['reads_all_tx']!s:5s} UMI lib {p['umi_library']:.0f} "
                  f"target {p['umi_target']} rest {p['umi_rest']} {fc} padj {pa} "
                  f"floor={p['passes_floor']} sig={p['sig']} (no floor {p['sig_without_floor']})")
    res["masked_probes"] = int(d.masked.sum())
    res["masked_genes"] = int(d.loc[d.masked, "gene_id"].nunique())
    nofloor = d.probe_sig_nofloor & ~d.gene_sig.fillna(False).astype(bool)
    res["masked_probes_without_floor"] = int(nofloor.sum())
    res["masked_genes_without_floor"] = int(d.loc[nofloor, "gene_id"].nunique())
    res["probe_de"] = int(d.probe_sig.sum())
    res["probe_de_without_floor"] = int(d.probe_sig_nofloor.sum())
    json.dump(res, open(o.out, "w"), indent=1)
    print(f"wrote {o.out}")


if __name__ == "__main__":
    main()
