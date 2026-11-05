#!/usr/bin/env python3
"""Figure 1: runtime and peak memory.

Runtime is the minimum over each series' timed reps; peak RSS is the maximum VmHWM
proc-tree sum (GiB). The values are read from the timing tables:
simpleaf from step 05 (results.tsv: dataset, wall_min_s, vmhwm_max_gib) and cyto / Cell Ranger
from step 05b (results.tsv: dataset, tool, wall_min_s,
vmhwm_max_gib). Every dataset must be present for every tool.

Style follows COMBINE-lab/alevin-fry-paper-scripts: seaborn `paper` context, `white`
style, despined, and that repo's palette.

Usage: make_fig1.py --simpleaf TSV --others TSV --out PNG --values-out TSV --axis-out JSON
  --axis-out  {"break_s": the runtime axis break} (the Figure 1 caption)
"""
import argparse
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ap = argparse.ArgumentParser()
ap.add_argument("--simpleaf", required=True); ap.add_argument("--others", required=True)
ap.add_argument("--out", required=True); ap.add_argument("--values-out", required=True)
ap.add_argument("--axis-out", required=True)
o = ap.parse_args()

PALETTE = ["#00AFBB", "#E7B800", "#D55E00"]
TOOLS = ["simpleaf", "cyto", "Cell Ranger"]
TOOL_C = dict(zip(TOOLS, PALETTE))
# display order and names of the five libraries (Table 1 of the manuscript)
LIBS = [("Flexv2_16k_Human", "16k human PBMC\nv2 4plex"),
        ("16k_mouse_Spleno", "16k mouse splenocyte\nv2 4plex"),
        ("40k_DTC", "40k human DTC\nv1 4plex"),
        ("16plex_human", "128k human PBMC\nv1 16plex"),
        ("320k_GEMX", "320k human PBMC\nv1 16plex")]
LIB = [n for _, n in LIBS]


def need(df, cols, what):
    miss = [c for c in cols if c not in df.columns]
    if miss:
        raise SystemExit(f"{what}: missing columns {miss} (have {list(df.columns)})")


sa = pd.read_csv(o.simpleaf, sep="\t"); need(sa, ["dataset", "wall_min_s", "vmhwm_max_gib"], o.simpleaf)
sa["tool"] = "simpleaf"
ot = pd.read_csv(o.others, sep="\t"); need(ot, ["dataset", "tool", "wall_min_s", "vmhwm_max_gib"], o.others)
ot["tool"] = ot["tool"].replace({"cellranger": "Cell Ranger"})
allv = pd.concat([sa[["dataset", "tool", "wall_min_s", "vmhwm_max_gib"]],
                  ot[["dataset", "tool", "wall_min_s", "vmhwm_max_gib"]]])
WALL, RSS = {}, {}
for tool in TOOLS:
    t = allv[allv.tool == tool].set_index("dataset")
    miss = [d for d, _ in LIBS if d not in t.index]
    if miss:
        raise SystemExit(f"{tool}: no timing for {miss}")
    WALL[tool] = [float(t.loc[d, "wall_min_s"]) for d, _ in LIBS]
    RSS[tool] = [float(t.loc[d, "vmhwm_max_gib"]) for d, _ in LIBS]
allv.to_csv(o.values_out, sep="\t", index=False)

# (value in axis units, label). Runtime lines are minutes expressed in seconds.
RSS_MARKS = [(8, "8 GiB"), (32, "32 GiB"), (128, "128 GiB")]

sns.set_context("paper")
sns.set_style("white")

# Runtime is linear with a break: simpleaf and cyto share the lower segment, Cell Ranger sits
# in the upper one. A log axis put all three within a factor of a few of each other and hid
# how far apart they are.
# The segments meet at the first whole 10 min above the slowest simpleaf / cyto run, so no bar
# crosses the boundary and no value is skipped
# -- only the scale changes. It must stay below Cell Ranger's fastest run.
LOWER_MAX = max(max(WALL["simpleaf"]), max(WALL["cyto"]))
BREAK = 600 * (math.floor(LOWER_MAX / 600) + 1)
assert LOWER_MAX < BREAK < min(WALL["Cell Ranger"]), f"no 10-min break between {LOWER_MAX} s and Cell Ranger's {min(WALL['Cell Ranger'])} s"
BREAK_LO = BREAK_HI = BREAK
TOP = 100 * math.ceil(max(WALL["Cell Ranger"]) / 100)

fig = plt.figure(figsize=(9.4, 4.2))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.25], hspace=0.22, wspace=0.34)
ax_hi = fig.add_subplot(gs[0, 0])
ax_lo = fig.add_subplot(gs[1, 0], sharex=ax_hi)
ax_rss = fig.add_subplot(gs[:, 1])
x = np.arange(len(LIB))
w = 0.26

for ax in (ax_hi, ax_lo):
    for k, tool in enumerate(TOOLS):
        v = np.array(WALL[tool], dtype=float)
        # in the upper segment, draw only the bars that reach it: a bar clipped to zero height
        # still paints its white edge along the boundary and reads as a gap
        if ax is ax_hi:
            v = np.where(v > BREAK, v, np.nan)
        ax.bar(x + (k - 1) * w, v, w, label=tool if ax is ax_hi else None,
               color=TOOL_C[tool], edgecolor="white", linewidth=0.5, zorder=3)
# both edges of the break are labeled, so the skipped range is readable
ax_lo.set_yticks([0, 1000, BREAK_LO])
ax_lo.set_ylim(0, BREAK_LO)
ax_hi.set_yticks([BREAK_HI, 10000, 20000])
ax_hi.set_ylim(BREAK_HI, TOP)

# minute marks, each drawn in whichever segment contains it
MINUTES = [1, 5, 10, 20, 30, 100, 200, 400]
for ax in (ax_lo, ax_hi):
    y0, y1 = ax.get_ylim()
    for mins in MINUTES:
        val = mins * 60
        if not (y0 < val <= y1):
            continue
        ax.axhline(val, color="0.55", linestyle=(0, (4, 3)), linewidth=0.7, zorder=1)
        ax.annotate("%d min" % mins, xy=(1.0, val), xycoords=("axes fraction", "data"),
                    xytext=(3, 0), textcoords="offset points", ha="left", va="center",
                    fontsize=7.5, color="0.4", zorder=4, annotation_clip=False)

# the break itself: the facing spines go (again at the end, after despine, which puts the
# bottom spine back), and the gap between the segments marks it
ax_hi.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
# the break marks, on the left spine only
d = 0.012
for ax, yy in ((ax_hi, 0), (ax_lo, 1)):
    ax.plot((-d, +d), (yy - d * 2.2, yy + d * 2.2), transform=ax.transAxes,
            color="0.3", lw=0.9, clip_on=False)

ax_lo.set_ylabel("time (s)")
ax_lo.yaxis.set_label_coords(-0.11, 1.0)
ax_hi.set_title("a  Runtime", loc="left", fontsize=10)
ax_lo.set_xticks(x)
ax_lo.set_xticklabels(LIB, rotation=45, ha="right")
# the legend sits above the upper segment, beside the panel title: inside it, it lands on
# the minute marks
ax_hi.legend(frameon=False, ncol=3, loc="lower right", bbox_to_anchor=(1.0, 1.0),
             fontsize=7.5, handlelength=1.1, columnspacing=1.1, handletextpad=0.5,
             borderpad=0.0, borderaxespad=0.2)

for k, tool in enumerate(TOOLS):
    ax_rss.bar(x + (k - 1) * w, RSS[tool], w, color=TOOL_C[tool], edgecolor="white",
               linewidth=0.5, zorder=3)
ax_rss.set_yscale("log")
ax_rss.set_ylim(bottom=1.4)
for val, lab in RSS_MARKS:
    ax_rss.axhline(val, color="0.55", linestyle=(0, (4, 3)), linewidth=0.7, zorder=1)
    ax_rss.annotate(lab, xy=(1.0, val), xycoords=("axes fraction", "data"), xytext=(3, 0),
                    textcoords="offset points", ha="left", va="center", fontsize=7.5,
                    color="0.4", zorder=4, annotation_clip=False)
ax_rss.set_ylabel("peak RSS (GiB)")
ax_rss.set_title("b  Peak memory", loc="left", fontsize=10)
ax_rss.set_xticks(x)
ax_rss.set_xticklabels(LIB, rotation=45, ha="right")

sns.despine(ax=ax_lo); sns.despine(ax=ax_hi); sns.despine(ax=ax_rss)
# thinner spines and ticks: at this bar width the default weight reads as a black band
for ax in (ax_lo, ax_hi, ax_rss):
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)
    ax.tick_params(width=0.6, length=3)
# the facing spines of the two runtime segments, hidden last so nothing restores them
ax_hi.spines["bottom"].set_visible(False)
ax_lo.spines["top"].set_visible(False)
fig.subplots_adjust(left=0.085, right=0.935, top=0.9, bottom=0.28)
print("break: lower 0-%d s, upper %d-%d s" % (BREAK_LO, BREAK_HI, TOP))
import json
json.dump({"break_s": BREAK}, open(o.axis_out, "w"))
fig.savefig(o.out, dpi=300)
print("wrote", o.out)
