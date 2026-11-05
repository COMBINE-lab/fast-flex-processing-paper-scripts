#!/usr/bin/env python3
"""Figure 2 layout, plotting only. Reads the cache written by umap_probe_plot.py and the
probe-relocation premrna table of the story library (panel c), so the layout can be changed
without re-reading the probe matrices.

The Figure 2 script (step 12 g).

Usage: fig2_plot.py --cache PKL --premrna TSV --out PNG

--cache is the step-12 output 12_probe_biology/figures/<ds>/fig_umap_cache.pkl and --premrna is
12_probe_biology/results/<ds>/premrna.tsv; step 12 calls this script with both.
"""
import argparse, pickle
import numpy as np, pandas as pd, seaborn as sns
import matplotlib; matplotlib.use("Agg")
import matplotlib.lines
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch, Polygon

from umap_orient import label_down

ap = argparse.ArgumentParser()
ap.add_argument("--cache", required=True); ap.add_argument("--premrna", required=True)
ap.add_argument("--out", required=True)
o = ap.parse_args()

C = pickle.load(open(o.cache, "rb"))
xy, lab, P4on, PD = C["xy"], C["lab"], C["P4on"], C["PD"]
GROUPS, H, NOTE, LABEL, GENE = C["GROUPS"], C["H"], C["NOTE"], C["label"], C["gene"]
# panel a: the two FIG2_PANEL_A examples (dilution left, cancelling right), each with its probes
# that pass the 10-UMI floor, shares per donor from the cache (umap_probe_plot.py)
PANEL_A = C["PANEL_A"]
GROUPS = [(e["gene"], e["label"], [p["probe"] for p in e["probes"] if p.get("floor", True)]) for e in PANEL_A]
PD = {(e["gene"], p["probe"]): list(zip(p["in"], p["out"])) for e in PANEL_A for p in e["probes"]}
# panel b: the outlined label points straight down (umap_orient.py, the same rotation
# Supp. Figure 2 uses); the hull outlines rotate with it
_rot = label_down(xy, lab, LABEL)
xy = _rot(xy)
H = {cl: (_rot(h), _rot(c)) for cl, (h, c) in H.items()}
# panel c: log(within-exon / junction-spanning probe signal) per label, as the former
# the junction/within-exon panel (probe_figures.py fig_premrna). Rule 2: the genes are those with a value in
# every label of >= MIN_LABEL_SAMPLES samples (pseudobulk units); a label with fewer samples is
# still drawn but does not choose the genes. The bar is the label's mean over those genes where it
# has a value, the dots are per-sample means over the same genes
MIN_LABEL_SAMPLES = 3
R = pd.read_csv(o.premrna, sep="\t")
_gl = R.groupby(["gene_id", "label"]).log_exon_over_junction.mean().unstack()
_nu = R.groupby("label")["sample"].nunique()
_g = _gl.loc[_gl[[c for c in _gl.columns if _nu[c] >= MIN_LABEL_SAMPLES]].notna().all(axis=1)]
_ps = (R[R.gene_id.isin(_g.index)].groupby(["label", "sample"]).log_exon_over_junction.mean()
       .unstack())
# shown as log(junction / within-exon): the negative of the table's log_exon_over_junction
CROWS = sorted(((c, -_g[c].dropna().mean(), list(-_ps.loc[c].dropna())) for c in _g.columns),
               key=lambda r: r[1])
N_PAIRED = _g.shape[0]   # the caption's gene count (rule 2)
SHORT = {"germinal_centre_B+proliferating": "germinal center B",
         "Tcm/Naive helper T cells": "naive helper T",
         "Tem/Trm cytotoxic T cells": "effector cytotoxic T",
         # the 16plex panel c has 15 labels: the long monocyte names ran into panel b's colorbar
         "Non-classical monocytes": "Non-classical mono.", "Classical monocytes": "Classical mono."}
SHORT = {c: SHORT.get(c, c.replace(" cells", "")) for c, _, _ in CROWS}
# the floor is light grey, not white: zero-count cells (a quarter of the embedding) must stay visible
CMAP = LinearSegmentedColormap.from_list(
    "af", ["#e6e6e6", "#b8e8ed", "#00AFBB", "#00727d", "#013238"])
AC = "#D55E00"; TEAL = "#00AFBB"

sns.set_context("paper"); sns.set_style("white")
# Explicit axes boxes, no tight_layout: panel a needs room below for the gene rules, b is a
# square embedding with a colorbar on its right, and c's cell-type labels need their own gutter.
fig = plt.figure(figsize=(7.2, 2.9))
_lo, _hi = xy.min(axis=0), xy.max(axis=0); _pad = 0.02 * (_hi - _lo)
XL, YL = (_lo[0] - _pad[0], _hi[0] + _pad[0]), (_lo[1] - _pad[1], _hi[1] + _pad[1])
_asp = 2.9 / 7.2 * (XL[1] - XL[0]) / (YL[1] - YL[0])   # width / height in figure fractions
BH = min(0.700, 0.250 / _asp); BW = BH * _asp           # equal aspect; width capped so the colorbar clears c
BY = 0.130 + (0.700 - BH) / 2
ax_b = fig.add_axes([0.345, BY, BW, BH])
ax_c = fig.add_axes([0.825, 0.190, 0.160, 0.700])   # full label names need the wider gutter

def heat(a_, v, title, vmax, cl, bar, mark=True, vmin=0):
    v = np.asarray(v, float); o = np.argsort(v)
    s = a_.scatter(*xy[o].T, s=0.5, c=v[o], cmap=CMAP, vmin=vmin, vmax=vmax, lw=0, rasterized=True)
    cax = fig.add_axes([0.345 + BW + 0.008, BY + 0.2 * BH, 0.010, 0.6 * BH])
    cb = fig.colorbar(s, cax=cax)
    cb.ax.tick_params(labelsize=6, length=2); cb.outline.set_visible(False)
    cb.set_label(bar, fontsize=6.2)
    if mark:
        hull, cen = H[cl]
        a_.add_patch(Polygon(hull, closed=True, fill=False, lw=0.8, ls=(0, (3, 2)), ec=AC,
                             zorder=6))
        a_.annotate(NOTE[cl], xy=(cen[0], hull[:, 1].min()), xytext=(0, -3),   # below: the label sits at the bottom after the rotation
                    textcoords="offset points", ha="center", va="top", fontsize=6.4,
                    color=AC, zorder=7)
    a_.set_xlim(*XL); a_.set_ylim(*YL); a_.set_aspect("equal")
    a_.set_xticks([]); a_.set_yticks([])
    for sp_ in a_.spines.values(): sp_.set_visible(False)

it = lambda g: "$\\it{%s}$" % g
heat(ax_b, P4on, None, np.percentile(P4on, 99), LABEL, "counts per 10k UMI")

a3 = ax_c
# drawn as the junction/within-exon panel (probe_figures.py fig_premrna): a bar per label at the
# gene mean, a dot per sample, a line at 0
for k, (cc, pooled, per) in enumerate(CROWS):
    a3.barh(k, pooled, color="#8fb8d8", height=0.62, edgecolor="#2f5b7c", linewidth=0.5, zorder=2)
    a3.scatter(per, np.full(len(per), k), s=2.5, color="#20303c", linewidth=0, zorder=3)
a3.axvline(0, color="#666", lw=0.6, zorder=1)
a3.set_yticks(np.arange(len(CROWS)))
a3.set_yticklabels([c.replace("Tcm/Naive", "Naive") for c, _, _ in CROWS], fontsize=5.2)
for t_ in a3.get_yticklabels():   # the naive cell types, the story of panel c, in bold
    if t_.get_text().startswith("Naive"): t_.set_fontweight("bold")
a3.set_ylim(-0.7, len(CROWS) - 0.3)
a3.set_xlabel("log(junction / within-exon)", fontsize=6.4)   # {N_PAIRED} genes (rule 2): caption
a3.tick_params(axis="x", labelsize=6); a3.tick_params(axis="y", length=0)
for sp_ in ("top", "right"):
    a3.spines[sp_].set_visible(False)

GROUPS = [(gene, cluster,
           sorted(probes, key=lambda k: np.mean([t[0] + t[1] for t in PD.get((gene, k), [(0, 0)])])))
          for gene, cluster, probes in GROUPS]
# one axis per gene inside panel a's box, widths proportional to the probe count
X0, W, GAP = 0.062, 0.235, 0.048
# DE probes (padj < 0.05, the 'sig' flag of the panel-a json) are bold and starred; the star is
# part of the tick text so the gene-rule extent below measures it
SIG = {(e["gene"], p["probe"]): p["sig"] for e in PANEL_A for p in e["probes"]}
AX_A = []
n = [len(p) for _, _, p in GROUPS]; unit = (W - GAP * (len(n) - 1)) / sum(n); x = X0
for gi, (gene, cluster, probes) in enumerate(GROUPS):
    a_ = fig.add_axes([x, 0.300, unit * n[gi], 0.590]); x += unit * n[gi] + GAP
    xt, xl = [], []
    for i, k in enumerate(probes):
        v = PD[(gene, k)]
        ins = np.array([t[0] for t in v]); out = np.array([t[1] for t in v])
        for xs, vv, cc_ in ((i - 0.16, out, "0.55"), (i + 0.16, ins, TEAL)):   # a box per side, 1.5 IQR whiskers
            a_.boxplot(vv, positions=[xs], widths=0.26, patch_artist=True, showfliers=True,
                       boxprops=dict(facecolor=cc_, edgecolor="0.25", lw=0.5),
                       medianprops=dict(color="0.15", lw=0.8), whiskerprops=dict(color="0.25", lw=0.5),
                       capprops=dict(color="0.25", lw=0.5),
                       flierprops=dict(marker="o", ms=1.5, mfc=cc_, mec="none"))
        xt.append(i); xl.append(k + ("*" if SIG[(gene, k)] else ""))
    a_.set_xticks(xt)
    tl = a_.set_xticklabels(xl, fontsize=5.2, rotation=55, ha="right", rotation_mode="anchor")
    for t_, k in zip(tl, probes):
        if SIG[(gene, k)]: t_.set_fontweight("bold")
    a_.set_xlim(-0.6, len(probes) - 0.4)
    a_.tick_params(axis="y", labelsize=6)
    AX_A.append((a_, gene))
    if gi == 0:
        a_.set_ylabel("share of the gene", fontsize=6.6)
    if gi == len(GROUPS) - 1:   # the legend sits in the wider right axis, upper left, clear of the boxes
        a_.legend(handles=[Patch(fc="0.55", ec="0.25", lw=0.5, label="other cells"),
                           Patch(fc=TEAL, ec="0.25", lw=0.5, label="the cluster")],
                  frameon=False, fontsize=6, loc="upper left", handlelength=1.0,
                  handletextpad=0.3, labelspacing=0.3, borderpad=0.1)
    sns.despine(ax=a_)

# all three titles on one baseline, each at its panel's left edge -- c's edge is its labels,
# not its axes, so the titles are placed in figure coordinates rather than per axes
TITLE_Y = 0.955
for x, t in ((0.062, "a  probes of two markers"),
             (0.345, "b  " + it(GENE) + " 202/203 probes"),
             (0.665, "c  junction / within-exon signal")):
    fig.text(x, TITLE_Y, t, ha="left", va="top", fontsize=8)

# gene rules under the probe labels: the x-extent is the rendered extent of the gene's rotated
# tick labels (left edge of the first to right edge of the last), so the line spans all of them
fig.canvas.draw()
_inv = fig.transFigure.inverted()
for a_, gene in AX_A:
    bb = [_inv.transform(t.get_window_extent()) for t in a_.get_xticklabels()]
    x0, x1 = min(b[0, 0] for b in bb), max(b[1, 0] for b in bb)
    yb = min(b[0, 1] for b in bb) - 2 / 72 / 2.9                # 2 pt below the lowest label
    fig.add_artist(matplotlib.lines.Line2D([x0, x1], [yb, yb], color="0.3", lw=0.7,
                                           solid_capstyle="butt"))
    fig.text((x0 + x1) / 2, yb - 3 / 72 / 2.9, gene, ha="center", va="top", fontsize=7,
             style="italic")
fig.savefig(o.out, dpi=400)
print(f"panel c: {N_PAIRED} genes (a value in every label with >= {MIN_LABEL_SAMPLES} samples), "
      f"{len(CROWS)} labels")
print("wrote", o.out)
