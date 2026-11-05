"""Orientation of the Figure 2 embedding, shared by Figure 2 panel b (fig2_plot.py) and
Supp. Figure 2 (fig_umap_labels.py), so the two figures always show the same orientation.

Rotation only, no reflection: the vector from the embedding centre (mean of all cells) to the
centroid of the `label` cells points straight down. Deterministic: depends only on xy and lab.
"""
import numpy as np


def label_down(xy, lab, label):
    """Return a function that maps points (n x 2) into the rotated frame."""
    xy = np.asarray(xy)
    m = np.asarray(lab) == label
    if not m.any():
        raise SystemExit(f"{label!r} not among the embedding labels")
    c0 = xy.mean(axis=0)
    v = xy[m].mean(axis=0) - c0
    t = -np.pi / 2 - np.arctan2(v[1], v[0])
    R = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    return lambda p: (np.asarray(p) - c0) @ R.T + c0
