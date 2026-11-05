#!/usr/bin/env python3
"""
Probe-vs-gene differential expression, both via PyDESeq2 (the sc-best-practices book's tool).

The claim this answers: aggregating probe counts to gene level HIDES differential expression
that is real at the probe level. A gene whose two probes move in opposite directions has a flat
total, so a gene-level test sees nothing while a probe-level test sees both probes move.

Method, identical at both granularities so the two are directly comparable:
  - Pseudobulk raw integer counts, paired target-vs-rest: for each donor, one "target" column
    (the cluster's cells) and one "rest" column (all that donor's other cells summed).
  - PyDESeq2 with design ``~ donor + grp`` -- donor as a blocking factor, grp = target/rest,
    exactly the book's multi-donor pseudobulk design (median-of-ratios size factors,
    per-feature negative-binomial dispersion, Wald test, BH-adjusted p).
  - Run once on the gene pseudobulk (probe matrix summed to genes) and once on the probe
    pseudobulk. A probe with padj < 0.05 whose gene has padj >= 0.05 is DE that gene-level
    aggregation masked.

Gene and probe DESeq2 runs each compute their own size factors and dispersions, as they would
if run standalone, so the discordance is a genuine disagreement between two proper analyses,
not an artefact of shared normalisation.

Minimum-UMI floor (Methods: "Probe testing was restricted to probes with at least 10 UMIs in both
groups"): only probes with at least --min-probe-umi UMI summed over the target columns AND summed
over the rest columns enter the probe-level DESeq2 run, so the size factors, dispersions and the
Benjamini-Hochberg correction cover the tested probes only. Untested probes stay in the table with
an empty padj and passes_floor False. Without the floor a probe with ~0 counts in the target is
"significant" through perfect separation (three PHACTR1 probes with 0-1 naive-B UMI were x4-x9).

Usage:
  pseudobulk_de.py --pb PB --target LABEL --out TSV --min-probe-umi N
  pseudobulk_de.py --selftest
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

MIN_TARGET_CELLS = 25      # a donor needs at least this many target cells to contribute a pair
MIN_DONORS = 3             # need at least this many paired donors to fit ~ donor + grp
ALPHA = 0.05


def _deseq(counts, meta, n_cpus=8):
    """counts: units x features int DataFrame; meta has 'donor' and 'grp'. -> results_df."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    import logging
    logging.getLogger("pydeseq2").setLevel(logging.ERROR)

    # drop features that are all-zero (DESeq2 cannot fit them, and they are not DE by anyone)
    nz = counts.columns[(counts.values.sum(axis=0) > 0)]
    counts = counts[nz]
    design = "~ donor + grp" if meta["donor"].nunique() > 1 else "~ grp"
    dds = DeseqDataSet(counts=counts, metadata=meta, design=design,
                       refit_cooks=True, quiet=True, n_cpus=n_cpus)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["grp", "target", "rest"], quiet=True, n_cpus=n_cpus)
    st.summary()
    return st.results_df


def build_pairs(pb, target):
    """-> (probe_counts units x probes int, meta) for the paired target-vs-rest design."""
    M = np.load(os.path.join(pb, "probe_counts.npy"))
    meta = pd.read_csv(os.path.join(pb, "units.tsv"), sep="\t")
    var = pd.read_csv(os.path.join(pb, "probes.tsv"), sep="\t", index_col=0)
    if target not in set(meta.label):
        raise SystemExit(f"{target!r} not in {sorted(set(meta.label))}")

    donors = sorted(meta["sample"].unique())
    rows, idx = [], []
    for d in donors:
        di = meta.index[meta["sample"] == d]
        tgt = di[meta.loc[di, "label"] == target]
        rest = di[meta.loc[di, "label"] != target]
        if len(tgt) == 0 or len(rest) == 0:
            continue
        if meta.loc[tgt, "n_cells"].sum() < MIN_TARGET_CELLS:
            continue
        rows.append(M[tgt].sum(axis=0)); idx.append(f"{d}__target")
        rows.append(M[rest].sum(axis=0)); idx.append(f"{d}__rest")
    if len({i.rsplit("__", 1)[0] for i in idx}) < MIN_DONORS:
        raise SystemExit(f"only {len(idx)//2} donors have >= {MIN_TARGET_CELLS} {target} cells")
    # sample names may themselves contain "__"; split on the LAST separator
    counts = pd.DataFrame(np.vstack(rows).astype(int), index=idx, columns=var.index)
    m = pd.DataFrame({"donor": [i.rsplit("__", 1)[0] for i in idx],
                      "grp": [i.rsplit("__", 1)[1] for i in idx]}, index=idx)
    return counts, m, var


def n_paired_donors(pb, target):
    """Donors build_pairs would use: the eligibility test the step applies before calling."""
    meta = pd.read_csv(os.path.join(pb, "units.tsv"), sep="\t")
    t = meta[meta.label == target].groupby("sample").n_cells.sum()
    has_rest = set(meta[meta.label != target]["sample"])
    return int(sum(1 for s_, n in t.items() if n >= MIN_TARGET_CELLS and s_ in has_rest))


def probe_deseq(counts, meta, min_probe_umi, n_cpus=8):
    """Probe-level DESeq2 on the probes with >= min_probe_umi UMI over the target columns AND over
    the rest columns only, so size factors, dispersions and the BH correction cover the tested
    probes. -> probe_-prefixed results for EVERY probe (untested probes: all NaN)."""
    tgt = meta.index[meta.grp == "target"]; rst = meta.index[meta.grp == "rest"]
    floor = (counts.loc[tgt].sum(axis=0) >= min_probe_umi) & (counts.loc[rst].sum(axis=0) >= min_probe_umi)
    return _deseq(counts.loc[:, floor.values], meta, n_cpus).add_prefix("probe_").reindex(counts.columns)


def probe_vs_gene(pb, target, out, min_probe_umi, n_cpus=8):
    counts, meta, var = build_pairs(pb, target)
    gene_of = var["gene_id"].reindex(counts.columns).values

    # gene pseudobulk = probe columns summed per gene, same units (transpose: groupby(axis=1)
    # is deprecated in pandas 2.x)
    gcounts = counts.T.groupby(gene_of).sum().T
    print(f"[de] {target}: {counts.shape[0]//2} paired donors, "
          f"{counts.shape[1]} probes / {gcounts.shape[1]} genes", flush=True)

    gres = _deseq(gcounts, meta, n_cpus).add_prefix("gene_")
    pres = probe_deseq(counts, meta, min_probe_umi, n_cpus)

    df = pres.join(pd.Series(gene_of, index=counts.columns, name="gene_id"))
    df = df.join(gres, on="gene_id")
    for c in ("gene_name", "genomic_class", "n_tx_hit", "n_tx_gene", "tx_set"):
        if c in var.columns:
            df[c] = var[c].reindex(df.index).values

    tgt = meta.index[meta.grp == "target"]; rst = meta.index[meta.grp == "rest"]
    df["umi_target"] = counts.loc[tgt].sum(axis=0).reindex(df.index).values
    df["umi_rest"] = counts.loc[rst].sum(axis=0).reindex(df.index).values
    df["passes_floor"] = (df["umi_target"] >= min_probe_umi) & (df["umi_rest"] >= min_probe_umi)
    df["probe_sig_nofloor"] = df["probe_padj"] < ALPHA
    df["probe_sig"] = df["probe_sig_nofloor"] & df["passes_floor"]
    df["gene_sig"] = df["gene_padj"] < ALPHA
    df["masked"] = df["probe_sig"] & (~df["gene_sig"].fillna(False))
    df = df.sort_values(["masked", "probe_padj"], ascending=[False, True])
    df.index.name = "probe_id"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    df.to_csv(out, sep="\t")

    n_masked = int(df["masked"].sum())
    n_masked_genes = df.loc[df["masked"], "gene_id"].nunique()
    print(f"[de] probe DE: {int(df.probe_sig.sum())} "
          f"({int(df.probe_sig_nofloor.sum())} before the {min_probe_umi}-UMI floor); "
          f"gene DE: {int(df.gene_sig.sum())}", flush=True)
    print(f"[de] MASKED (probe padj<{ALPHA}, gene padj>={ALPHA}): "
          f"{n_masked} probes on {n_masked_genes} genes -> {out}", flush=True)
    show = df[df["masked"]].head(15)
    cols = [c for c in ("gene_name", "genomic_class", "n_tx_hit", "n_tx_gene",
                        "probe_log2FoldChange", "probe_padj", "gene_log2FoldChange",
                        "gene_padj") if c in show.columns]
    if len(show):
        with pd.option_context("display.width", 200, "display.max_columns", 20):
            print(show[cols].to_string(), flush=True)
    return df


def selftest():
    """Canonical masking case: one gene, two probes moving oppositely, gene flat."""
    rng = np.random.default_rng(0)
    ndon = 6
    idx, donor, grp, rows = [], [], [], []
    # features: G_up, G_down (one gene G), then 40 null probes each its own gene
    feat = ["G|GENE|up", "G|GENE|down"] + [f"N{i}|N{i}|x" for i in range(40)] + ["L|LOW|x"]
    for d in range(ndon):
        for g in ("target", "rest"):
            base = rng.poisson(300, len(feat)).astype(float)
            if g == "target":
                # additive opposite shifts so the GENE TOTAL is unchanged (a x4/-4 pair does
                # not cancel in counts: 300*4 + 300/4 = 1275, not 600). up+down must stay 2*300.
                base[0] = base[0] + 200   # probe up   ~500
                base[1] = base[1] - 200   # probe down ~100  -> gene total ~ 600 either way
            base[-1] = 0.05 if g == "target" else 300   # a probe under the floor in the target
            idx.append(f"D{d}__{g}"); donor.append(f"D{d}"); grp.append(g)
            rows.append(rng.poisson(base))
    counts = pd.DataFrame(np.vstack(rows).astype(int), index=idx, columns=feat)
    meta = pd.DataFrame({"donor": donor, "grp": grp}, index=idx)
    gene_of = np.array([f.split("|")[0] for f in feat])
    gcounts = counts.T.groupby(gene_of).sum().T

    gres = _deseq(gcounts, meta, 4)
    pres = probe_deseq(counts, meta, 10, 4).rename(columns=lambda c: c[len("probe_"):])
    low_untested = bool(pres.loc["L|LOW|x"].isna().all())
    print(f"[selftest] probe under the 10-UMI floor untested (padj empty): {low_untested}")
    g_padj = gres.loc["G", "padj"]
    up_padj, dn_padj = pres.loc["G|GENE|up", "padj"], pres.loc["G|GENE|down", "padj"]
    print(f"[selftest] gene G padj={g_padj:.3g} (expect NOT sig)")
    print(f"[selftest] probe up padj={up_padj:.3g}, down padj={dn_padj:.3g} (expect BOTH sig)")
    null_sig = int((pres.loc[[f for f in feat if f.startswith('N')], 'padj'] < ALPHA).sum())
    print(f"[selftest] null probes called DE: {null_sig}/40 (expect ~0)")
    ok = (g_padj >= ALPHA) and (up_padj < ALPHA) and (dn_padj < ALPHA) and null_sig <= 2 and low_untested
    print(f"[selftest] {'PASS' if ok else 'FAIL'}: gene masked, both probes recovered, floor probe untested")
    sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb"); ap.add_argument("--target"); ap.add_argument("--out")
    ap.add_argument("--n-cpus", type=int, default=8)
    ap.add_argument("--min-probe-umi", type=int)
    ap.add_argument("--count-donors", action="store_true",
                    help="print the number of donors eligible for the paired design and exit")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
    if a.count_donors:
        print(n_paired_donors(a.pb, a.target), MIN_DONORS)
        return
    if not (a.pb and a.target and a.out and a.min_probe_umi is not None):
        ap.error("need --pb, --target, --out, --min-probe-umi (or --selftest)")
    probe_vs_gene(a.pb, a.target, a.out, a.min_probe_umi, a.n_cpus)


if __name__ == "__main__":
    main()
