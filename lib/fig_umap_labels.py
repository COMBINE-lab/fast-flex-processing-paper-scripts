#!/usr/bin/env python3
"""Supp. Figure 2: cell labels of the Figure 2 embedding, with the number of cells per label.

Every label with at least --min-cells cells is colored, the rest are greyed, and the
--highlight label (the one Figure 2 panel b outlines) is drawn on top.

Usage: fig_umap_labels.py --embedding H5AD --out PNG --highlight "Naive B cells" [--min-cells 500]
"""
import argparse

import anndata as ad
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from umap_orient import label_down

ap = argparse.ArgumentParser()
ap.add_argument("--embedding", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--highlight", required=True); ap.add_argument("--min-cells", type=int, default=500)
o = ap.parse_args()

G = ad.read_h5ad(o.embedding)
xy = G.obsm["X_umap"]; lab = G.obs["label"].astype(str).values
xy = label_down(xy, lab, o.highlight)(xy)   # same orientation as Figure 2 panel b
labs, cnt = np.unique(lab, return_counts=True)
order = [l for _, l in sorted(zip(-cnt, labs))]
counts = dict(zip(labs, cnt))
named = [l for l in order if counts[l] >= o.min_cells]
if o.highlight not in counts:
    raise SystemExit(f"{o.highlight!r} not among the embedding labels")
pal = dict(zip(named, sns.color_palette("husl", len(named))))

sns.set_context("paper"); sns.set_style("white")
fig, ax = plt.subplots(figsize=(5.4, 4.2))
other = ~np.isin(lab, named)
if other.any():
    ax.scatter(*xy[other].T, s=0.4, c="0.88", lw=0, rasterized=True,
               label=f"other labels ({other.sum():,})")
for l in named:
    m = lab == l
    ax.scatter(*xy[m].T, s=0.4, c=[pal[l]], lw=0, rasterized=True, label=f"{l} ({m.sum():,})",
               zorder=6 if l == o.highlight else 2)
ax.legend(frameon=False, fontsize=5.4, loc="center left", bbox_to_anchor=(1.0, 0.5),
          markerscale=4, handletextpad=0.2, labelspacing=0.3)
ax.set_title(f"{G.n_obs:,} cells, {len(G.uns['donors'])} donors", loc="left", fontsize=8)
ax.set_xticks([]); ax.set_yticks([])
for sp_ in ax.spines.values():
    sp_.set_visible(False)
fig.tight_layout(); fig.savefig(o.out, dpi=400)
print("wrote", o.out)
