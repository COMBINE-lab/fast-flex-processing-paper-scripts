#!/usr/bin/env python3
"""Figure 2 panel a, the simple way: per cell type, take its CellTypist markers that are flat at
gene level and have a differentially expressed probe, and classify each as cancelling or diluted.

Input: the gene-vs-probe DE tables of step 12c (one TSV per label, columns as written by
lib/pseudobulk_de.py), one directory per dataset.

  cancelling: a significant probe, and another floor-passing probe moving the other way by >= --min-opp-lfc.
  diluted:    >=1 significant probe that is not the gene's largest probe, while the largest
              probe is not significant and moves less than half as much (|log2fc|).
  flat gene:  gene not DE, gene_padj >= --flat-padj (default 0.05, the paper's masked rule).
  marker:     CellTypist Immune_All_Low coefficient rank of the gene for the label's class <= --max-rank.

Usage: fig2_simple_examples.py --model PKL --out TSV --de DATASET=DIR [--de ...]
"""
import argparse, glob, os, re
import numpy as np, pandas as pd

LABEL_MAP = {  # labels that are not CellTypist classes themselves
    "proliferating": ["Cycling B cells", "Cycling DCs", "Cycling NK cells", "Cycling T cells", "Cycling monocytes"],
    "plasmablast": ["Plasmablasts"],
    "germinal_centre_B+proliferating": ["Germinal center B cells"],
    "Monocytes": ["Classical monocytes", "Non-classical monocytes"],
}


def ranks(model, label):
    cls = list(model.classifier.classes_)
    cs = [c for c in LABEL_MAP.get(label, [label]) if c in cls]
    if not cs:
        return None
    r = None
    for c in cs:
        order = np.argsort(-model.classifier.coef_[cls.index(c)])
        s = pd.Series(np.arange(1, len(order) + 1), index=np.asarray(model.features)[order])
        r = s if r is None else np.minimum(r, s.reindex(r.index))
    return r


def classify(g, alpha, min_opp):
    p = g[g.passes_floor.astype(bool)]
    if len(p) < 2:
        return None
    sig = p[p.probe_padj < alpha]
    if sig.empty:
        return None
    # cancelling (the paper's ZBTB16): a significant probe, and another floor-passing probe moving the
    # other way by a visible amount (it need not be significant itself)
    a = sig.loc[sig.probe_log2FoldChange.abs().idxmax()]
    opp = p[np.sign(p.probe_log2FoldChange) == -np.sign(a.probe_log2FoldChange)]
    opp = opp[opp.probe_log2FoldChange.abs() >= min_opp]
    if not opp.empty:
        b = opp.loc[opp.probe_log2FoldChange.abs().idxmax()]
        return "cancelling", a, b
    big = p.loc[p.probe_baseMean.idxmax()]
    movers = sig.drop(big.name, errors="ignore")
    if (not movers.empty and not (big.probe_padj < alpha)):
        m = movers.loc[movers.probe_log2FoldChange.abs().idxmax()]
        if abs(big.probe_log2FoldChange) <= 0.5 * abs(m.probe_log2FoldChange):
            return "diluted", m, big
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--de", nargs="+", required=True)
    ap.add_argument("--flat-padj", type=float, default=0.05); ap.add_argument("--flat-lfc", type=float, default=1.0)
    ap.add_argument("--min-opp-lfc", type=float, default=0.1)   # |log2fc| of the opposing probe (0.1 = 7%)
    ap.add_argument("--alpha", type=float, default=0.05); ap.add_argument("--max-rank", type=int, default=200)
    o = ap.parse_args()
    import celltypist
    model = celltypist.models.Model.load(o.model)
    rows = []
    for spec in o.de:
        ds, d = spec.split("=", 1)
        for f in sorted(glob.glob(os.path.join(d, "*.tsv"))):
            t = pd.read_csv(f, sep="\t")
            if t.empty or "gene_padj" not in t:
                continue
            label = t.attrs.get("label") or os.path.basename(f)[:-4].replace("_", " ")
            label = {"CD16+ NK cells": "CD16+ NK cells", "Tcm Naive cytotoxic T cells": "Tcm/Naive cytotoxic T cells",
                     "Tcm Naive helper T cells": "Tcm/Naive helper T cells", "Tem Temra cytotoxic T cells": "Tem/Temra cytotoxic T cells",
                     "Tem Effector helper T cells": "Tem/Effector helper T cells",
                     "germinal centre B+proliferating": "germinal_centre_B+proliferating"}.get(label, label)
            rk = ranks(model, label)
            # mouse: model features are human symbols; upper-case is the usual 1:1 approximation
            key = (lambda s: s.upper()) if "mouse" in ds else (lambda s: s)
            t = t[t.gene_tested.astype(bool)] if "gene_tested" in t else t
            flat = t[(t.gene_padj >= o.flat_padj) & (t.gene_log2FoldChange.abs() <= o.flat_lfc)]
            for gene, g in flat.groupby("gene_name"):
                r = rk.get(key(gene), np.nan) if rk is not None else np.nan
                if not (r <= o.max_rank):
                    continue
                c = classify(g, o.alpha, o.min_opp_lfc)
                if c is None:
                    continue
                kind, a, b = c
                rows.append(dict(dataset=ds, label=label, gene=gene, story=kind, marker_rank=int(r),
                                 gene_fc=2 ** g.gene_log2FoldChange.iloc[0], gene_padj=g.gene_padj.iloc[0],
                                 probe_a=a.probe_id.split("|")[-1], fc_a=2 ** a.probe_log2FoldChange, padj_a=a.probe_padj,
                                 umi_a=f"{int(a.umi_target)}/{int(a.umi_rest)}", tx_a=a.tx_set,
                                 probe_b=b.probe_id.split("|")[-1], fc_b=2 ** b.probe_log2FoldChange, padj_b=b.probe_padj,
                                 umi_b=f"{int(b.umi_target)}/{int(b.umi_rest)}", tx_b=b.tx_set,
                                 n_probes=len(g), n_tx_gene=int(a.n_tx_gene)))
    out = pd.DataFrame(rows).sort_values(["story", "marker_rank", "gene_padj"], ascending=[True, True, False])
    out.to_csv(o.out, sep="\t", index=False)
    for k in ("cancelling", "diluted"):
        s = out[out.story == k]
        print(f"\n== {k}: {len(s)} flat markers (rank <= {o.max_rank}) with a DE probe")
        for _, r in s.head(12).iterrows():
            print(f"  {r.dataset:16s} {r.label:28s} {r.gene:10s} rank {r.marker_rank:4d}  gene x{r.gene_fc:.2f} (padj {r.gene_padj:.2f})"
                  f"  {r.probe_a} x{r.fc_a:.2f} (padj {r.padj_a:.1e}, {r.umi_a})  vs  {r.probe_b} x{r.fc_b:.2f} (padj {r.padj_b:.1e}, {r.umi_b})")


if __name__ == "__main__":
    main()
