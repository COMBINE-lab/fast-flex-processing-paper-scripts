#!/usr/bin/env python3
"""Supp. Figure 1: the pooled 128k human PBMC embedding, colored by each tool's cell-type labels.

All three panels draw the SAME cells at the SAME positions: Cell Ranger's pooled UMAP (step 10,
book variant). Panel a colors each cell by its simpleaf label, b by its cyto label, c by its
Cell Ranger label (supp1_tool_labels.py, step 11's labelling on each tool's own counts). A label
has the same color in every panel, so a cell whose color changes between panels is a cell the
tools label differently.

Usage: supp_fig1_pooled_umap.py --embedding OBS.TSV.GZ --labels TSV.GZ --out PNG [--values-out JSON] [--min-cells 500]
  --embedding  step-10 <ds>_cellranger_book.obs.tsv.gz (index sample_barcode; sample, umap1, umap2)
  --labels     one table, columns cell, tool, label (cell = sample_barcode)

--embedding is the step-10 output <ds>_cellranger_book.obs.tsv.gz; --labels comes from
reproduce/lib/supp1_tool_labels.py on the step-09 objects. Step 14 calls this script with both.
"""
import argparse

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TOOLS = [("simpleaf", "a  simpleaf"), ("cyto", "b  cyto"), ("cellranger", "c  Cell Ranger")]

ap = argparse.ArgumentParser()
ap.add_argument("--embedding", required=True); ap.add_argument("--labels", required=True)
ap.add_argument("--out", required=True); ap.add_argument("--min-cells", type=int, default=500)
ap.add_argument("--values-out", help="JSON: samples, cells and the share of cells labelled as by Cell Ranger, per tool")
o = ap.parse_args()

E = pd.read_csv(o.embedding, sep="\t", index_col=0)
L = pd.read_csv(o.labels, sep="\t").pivot(index="cell", columns="tool", values="label")
cells = E.index.intersection(L.dropna().index)
if len(cells) < 0.99 * len(E):
    raise SystemExit(f"only {len(cells)} of {len(E)} embedded cells carry a label in all three tools")
E, L = E.loc[cells], L.loc[cells]
xy = E[["umap1", "umap2"]].to_numpy()

# one palette over the labels any tool gives to >= min-cells cells; the rest are grey
cnt = pd.concat([L[t] for t, _ in TOOLS]).value_counts()
named = [l for l in cnt.index if cnt[l] / 3 >= o.min_cells]
pal = dict(zip(named, sns.color_palette("husl", len(named))))

sns.set_context("paper"); sns.set_style("white")
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.75))
lo, hi = xy.min(axis=0), xy.max(axis=0); pad = 0.02 * (hi - lo)
for ax, (tool, title) in zip(axs, TOOLS):
    lab = L[tool].values
    other = ~np.isin(lab, named)
    ax.scatter(*xy[other].T, s=0.15, c="0.88", lw=0, rasterized=True)
    for l in named:
        m = lab == l
        ax.scatter(*xy[m].T, s=0.15, c=[pal[l]], lw=0, rasterized=True)
    ax.set_title(title, loc="left", fontsize=8)
    ax.set_xlim(lo[0] - pad[0], hi[0] + pad[0]); ax.set_ylim(lo[1] - pad[1], hi[1] + pad[1])
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
handles = [plt.Line2D([], [], ls="", marker="o", ms=3.5, color=pal[l]) for l in named]
fig.legend(handles, [l.replace("Tcm/Naive", "Naive") for l in named], loc="lower center",
           ncol=5, frameon=False, fontsize=5.6, handletextpad=0.2, columnspacing=0.9)
fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.2, wspace=0.03)
fig.savefig(o.out, dpi=400)
same = {t: float((L[t] == L["cellranger"]).mean()) for t, _ in TOOLS}
print(f"wrote {o.out}: {len(cells):,} cells; same label as Cell Ranger: "
      + ", ".join(f"{t} {v:.3f}" for t, v in same.items()))
if o.values_out:
    import json
    json.dump({"n_samples": int(E["sample"].nunique()), "n_cells": int(len(cells)), "same_label_as_cellranger": same}, open(o.values_out, "w"), indent=1)
