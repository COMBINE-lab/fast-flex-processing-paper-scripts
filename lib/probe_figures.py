#!/usr/bin/env python3
"""
Figures of the probe-level analysis used by the manuscript.

  fig1_premrna_by_celltype.png   instrument B's log(exon/junction) per label,
                                 over the genes measured in every label, per-sample dots
  fig4_de_gene_verdicts.png      what the probe level says about each cluster's
                                 DE genes (from genede_join__*.tsv)

The results directory is an argument. verdict_summary.tsv is written next to the figure.

Usage: probe_figures.py --dataset DS --results DIR --outdir DIR
"""
import argparse, glob, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"figure.dpi": 130, "font.size": 8, "axes.spines.top": False,
                     "axes.spines.right": False})


def fig_premrna(ds, results, outdir):
    d = pd.read_csv(f"{results}/premrna.tsv", sep="\t")
    g = d.groupby(["gene_id", "label"]).log_exon_over_junction.mean().unstack().dropna()
    order = g.mean().sort_values().index
    # The per-sample dots must be computed on the SAME paired gene subset as the bars,
    # otherwise each label's dots average a different gene set from its bar and land off it.
    paired = d[d.gene_id.isin(g.index)]
    per_sample = (paired.groupby(["label", "sample"]).log_exon_over_junction.mean()
                  .unstack().reindex(order))

    fig, ax = plt.subplots(figsize=(6.4, 0.34 * len(order) + 1.3))
    y = np.arange(len(order))
    ax.barh(y, g.mean().reindex(order).values, color="#8fb8d8", height=0.62,
            edgecolor="#2f5b7c", linewidth=0.6, zorder=2)
    for k, lb in enumerate(order):
        vals = per_sample.loc[lb].dropna().values
        ax.scatter(vals, np.full(len(vals), k), s=9, color="#20303c", zorder=3,
                   linewidth=0)
    ax.axvline(0, color="#666", linewidth=0.7, zorder=1)
    ax.set_yticks(y); ax.set_yticklabels(order)
    ax.set_xlabel("log(within-exon / junction probe signal)   ->  higher = more pre-mRNA")
    ax.set_title(f"{ds}: splicing state by cell type\n"
                 f"isoform-set-matched probe pairs, {g.shape[0]} genes; dots = per-sample means",
                 fontsize=8.5, loc="left")
    fig.tight_layout(); fig.savefig(f"{outdir}/fig1_premrna_by_celltype.png"); plt.close(fig)
    print("  fig1 splicing state")


def fig_verdicts(ds, results, outdir):
    rows = []
    files = sorted(glob.glob(f"{results}/genede_join__*.tsv"))
    if not files:
        raise SystemExit(f"no genede_join__*.tsv under {results}")
    for f in files:
        lb = os.path.basename(f)[len("genede_join__"):-4].rstrip("_")
        j = pd.read_csv(f, sep="\t")
        sig = j[j.q < 0.05]
        if not len(sig):
            continue
        v = sig.probe_verdict.value_counts()
        rows.append({"label": lb, "de_genes": len(sig),
                     "concordant": v.get("concordant_expression", 0),
                     "isoform_switch": v.get("isoform_switch", 0),
                     "splicing_state": v.get("splicing_state", 0),
                     "switch_not_de": int(((j.q >= 0.05) &
                                           (j.probe_verdict == "isoform_switch")).sum())})
    if not rows:
        raise SystemExit(f"{ds}: no cluster has a DE gene; nothing to plot")
    t = pd.DataFrame(rows)
    t["total"] = t.isoform_switch + t.splicing_state + t.switch_not_de
    t = t[t.total > 0].sort_values("total")
    if not len(t):
        raise SystemExit(f"{ds}: no cluster has a probe-level verdict; nothing to plot")
    # Counts span three orders of magnitude across clusters (one cluster can carry >1000 while
    # most carry single digits), so a linear axis renders every cluster but the largest as a
    # zero-width sliver. symlog keeps the small clusters legible without hiding the big one.
    fig, ax = plt.subplots(figsize=(7.2, 0.36 * len(t) + 1.7))
    y = np.arange(len(t))
    ax.barh(y, t.isoform_switch, color="#c0603f", label="isoform switch (DE)")
    ax.barh(y, t.splicing_state, left=t.isoform_switch, color="#3f7ec0",
            label="splicing-state shift (DE)")
    ax.barh(y, t.switch_not_de, left=t.isoform_switch + t.splicing_state,
            color="#e0b04a", label="isoform switch, NOT DE")
    ax.set_xscale("symlog", linthresh=1)
    for k, r in enumerate(t.itertuples()):
        ax.text(r.total * 1.25 + 0.6, k, f"of {r.de_genes:,} DE", va="center",
                fontsize=6.5, color="#444")
    ax.set_yticks(y); ax.set_yticklabels(t.label, fontsize=7)
    ax.set_xlim(0, t.total.max() * 6)
    ax.set_xlabel("genes (symlog scale)")
    ax.legend(fontsize=6.5, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.16 / (0.36 * len(t) + 1.7) * 6 + 0.02), ncol=3)
    ax.set_title(f"{ds}: what the probe level adds to each cluster's DE genes\n"
                 f"gold = isoform switch with no gene-level DE, invisible to a gene pipeline",
                 fontsize=8.5, loc="left")
    fig.tight_layout(); fig.savefig(f"{outdir}/fig4_de_gene_verdicts.png",
                                    bbox_inches="tight"); plt.close(fig)
    print("  fig4 DE gene verdicts")
    t.to_csv(f"{outdir}/verdict_summary.tsv", sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--which", nargs="+", default=["premrna", "verdicts"],
                    choices=["premrna", "verdicts"])
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    print(f"[figures] {a.dataset}: {a.results} -> {a.outdir}")
    if "premrna" in a.which:
        fig_premrna(a.dataset, a.results, a.outdir)
    if "verdicts" in a.which:
        fig_verdicts(a.dataset, a.results, a.outdir)


if __name__ == "__main__":
    main()
