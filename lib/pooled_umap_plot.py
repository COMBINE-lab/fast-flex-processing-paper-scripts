#!/usr/bin/env python3
"""1 x 3 pooled-UMAP panels (simpleaf, cyto, Cell Ranger) for one pooled_umap.py variant.

Style follows flex-paper tools/umap_probe_plot.py: points s=0.5 rasterized, 6-8 pt text, left-aligned
panel titles ('a  simpleaf'), no axes, 400 dpi. Colors: Cell Ranger cluster c gets palette color c; every
simpleaf / cyto cluster takes the color of the Cell Ranger cluster it was Hungarian-matched to, so equal
colors are matched clusters; clusters left unmatched (u*) are drawn in grey.

Usage: pooled_umap_plot.py --dir DIR --dataset DS --variant basic|book --out PNG
"""
import argparse, json

import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scanpy.plotting.palettes import godsnot_102

NAMES = [("simpleaf", "simpleaf"), ("cyto", "cyto"), ("cellranger", "Cell Ranger")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True); ap.add_argument("--dataset", required=True)
    ap.add_argument("--variant", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    obs = {t: pd.read_csv(f"{a.dir}/{a.dataset}_{t}_{a.variant}.obs.tsv.gz", sep="\t", index_col=0,
                          dtype={"leiden": str, "matched": str}) for t, _ in NAMES}
    conc = json.load(open(f"{a.dir}/{a.dataset}_{a.variant}_concordance.json"))
    cr = obs["cellranger"]["leiden"]
    order = cr.value_counts().index.tolist()                 # largest Cell Ranger cluster -> first color
    color = {c: godsnot_102[i % len(godsnot_102)] for i, c in enumerate(order)}

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.55))
    for i, (ax, (t, title)) in enumerate(zip(axes, NAMES)):
        o = obs[t]
        xy = o[["umap1", "umap2"]].to_numpy()
        cols = np.array([color.get(m, "#bdbdbd") for m in o["matched"]])
        perm = np.random.default_rng(0).permutation(len(o))   # no cluster systematically drawn on top
        ax.scatter(xy[perm, 0], xy[perm, 1], s=0.5, c=cols[perm], lw=0, rasterized=True)
        for m, g in o.groupby("matched"):
            if m in color and len(g) >= 0.004 * len(o):
                med = g[["umap1", "umap2"]].median()
                ax.text(med.iloc[0], med.iloc[1], m, fontsize=6, ha="center", va="center", color="0.1")
        ax.set_title(f"{'abc'[i]}  {title}", loc="left", fontsize=8)
        sub = f"{o['leiden'].nunique()} clusters"
        if t != "cellranger":
            p = conc["pairs"][f"cellranger_vs_{t}"]
            sub += f"; vs Cell Ranger ARI {p['ARI']:.3f}, NMI {p['NMI']:.3f}"
        ax.text(0.0, -0.02, sub, transform=ax.transAxes, fontsize=6, ha="left", va="top")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    fig.text(0.995, 0.005, f"{conc['n_cells']:,} cells, {a.variant} pipeline, batch correction: "
             f"{conc['batch_correction'].split(' (')[0]}", fontsize=6, ha="right", va="bottom", color="0.35")
    fig.tight_layout(w_pad=1.0)
    fig.savefig(a.out, dpi=400)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
