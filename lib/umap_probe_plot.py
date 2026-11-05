#!/usr/bin/env python3
"""Figure 2 cache: every number panels a and b need, so the layout (fig2_plot.py) never re-reads
the probe matrices.

  a  for each --panel-a example "DATASET|CELL TYPE|GENE" (config FIG2_PANEL_A_DILUTED and
     FIG2_PANEL_A_CANCELLING; the two may come from different libraries and cell types): each
     probe's share of the gene, per donor of that library, in the cell type and in the donor's
     other cells (descriptive). Donors of the embedding's library follow --samples. Each example
     keeps its own donors: those with >= --min-cluster-cells cells of the type and >= --min-gene-umi
     gene UMI in the type and in the other cells (a subset of the embedding's donors).
  b  per cell of the embedding, the counts of the --gene probes whose isoform set contains every
     --tx transcript, pooled with the total UMI over the --k-pool nearest cells in PCA space,
     as counts per 10k UMI.

Paths and the story constants are arguments. Probes are selected by gene name, asserting their
number. Cells are labelled by the `label` column, the same as the pseudobulk DE. Panel c is read
from the allgene_ratio table by fig2_plot.py, not cached here.

Usage:
  umap_probe_plot.py --probe-root DIR --label-root DIR --dataset DS --embedding H5AD --out PKL
     --label "Naive B cells" --panel-a "320k_GEMX|Naive B cells|YBX3" "320k_GEMX|Naive B cells|ZBTB16"
     --gene PHACTR1
     --tx ENST00000379329 ENST00000379335 [--k-pool 100] [--hull-labels "Naive B cells" ...]
"""
import argparse
import os
import pickle
import sys

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from probe_instruments import probe_paths, read_labels


def island(xy, lab, cl, trim=70, pad=0.10):
    from scipy.spatial import ConvexHull
    p = xy[lab == cl]
    if len(p) < 3:
        raise SystemExit(f"label {cl!r} has {len(p)} cells in the embedding; cannot outline it")
    dd = np.linalg.norm(p - np.median(p, axis=0), axis=1)
    keep = p[dd <= np.percentile(dd, trim)]
    cen = keep.mean(axis=0)
    h = keep[ConvexHull(keep).vertices]
    return cen + (h - cen) * (1 + pad), cen


def load_panel_a(o):
    """-> PD {(gene, probe): [(share in cell type, share in other cells) per donor]}, probes per
    gene, donors per example, for every --panel-a example from its own library. A donor that fails
    a per-donor rule is left out and counted; an example with no donor left fails, naming the rule."""
    import anndata as ad
    PD, probes_of, donors = {}, {}, {}
    for ex in o.panel_a:
        ds, label, gene = ex.split("|")
        donors[ex] = []
        paths = probe_paths(o.probe_root, ds, o.samples if ds == o.dataset else None)
        few_cells, low_umi = [], []    # donors that fail the cell rule / the gene-UMI rule
        for d, path in paths:
            a = ad.read_h5ad(path)
            lb = read_labels(o.label_root, ds, d, o.drop_doublets)
            keep = [b for b in a.obs_names if b in lb.index]
            a = a[keep]; labs = lb.loc[keep].values
            gn = a.var["gene_name"].astype(str).values
            sel = np.flatnonzero(gn == gene)
            if not len(sel):
                raise SystemExit(f"{ex}: {gene} has no probes in the panel")
            probes_of[gene] = [s.split("|")[-1] for s in np.array(a.var_names)[sel]]
            m = labs == label
            if m.sum() < o.min_cluster_cells:
                few_cells.append(d)
                continue
            X = a.X.tocsr()
            Xin, Xout = X[m][:, sel], X[~m][:, sel]
            din, dout = Xin.sum(), Xout.sum()
            if din < o.min_gene_umi or dout < o.min_gene_umi:
                low_umi.append(d)
                continue
            donors[ex].append(d)
            pin = np.asarray(Xin.sum(axis=0)).ravel() / din
            pout = np.asarray(Xout.sum(axis=0)).ravel() / dout
            for k, pid in enumerate(probes_of[gene]):
                PD.setdefault((gene, pid), []).append((float(pin[k]), float(pout[k])))
        rules = (f"{len(few_cells)} with < {o.min_cluster_cells} cells of the type, {len(low_umi)} "
                 f"with < {o.min_gene_umi} {gene} UMI in the type or in the other cells")
        if not donors[ex]:
            raise SystemExit(f"{ex}: no donor of {len(paths)} passes the per-donor rules ({rules})")
        print(f"[panel a] {ex}: {len(donors[ex])} of {len(paths)} donors (left out: {rules})", flush=True)
    return PD, probes_of, donors


def check_donors(ex_donors, dataset, emb_donors):
    """FIG_UMAP_DONORS: panels a and b draw on one donor set. An example of the embedding's library
    must not use a donor outside the embedding (e.g. an embedding built for another donor set); it
    may use fewer, because load_panel_a's per-donor rules leave donors out per example."""
    for ex, dn in ex_donors.items():
        if ex.split("|")[0] == dataset and not set(dn) <= emb_donors:
            raise SystemExit(f"panel a {ex}: donors {sorted(set(dn) - emb_donors)} are not in the "
                             f"embedding (donors {sorted(emb_donors)})")


def main():
    import anndata as ad
    from sklearn.neighbors import NearestNeighbors
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe-root", required=True); ap.add_argument("--label-root", required=True)
    ap.add_argument("--dataset", required=True); ap.add_argument("--embedding", required=True)
    ap.add_argument("--samples", nargs="+", help="restrict panel a to these donors (tests)")
    ap.add_argument("--drop-doublets", action="store_true")
    ap.add_argument("--out", required=True); ap.add_argument("--label", required=True)
    ap.add_argument("--panel-a", nargs="+", required=True, metavar="DS|LABEL|GENE")
    ap.add_argument("--panel-a-json", nargs="+", required=True, metavar="ROLE=JSON",
                    help="story_numbers.py output per example (probe fold change, call, isoforms); "
                         "roles in --panel-a order")
    ap.add_argument("--gene", required=True); ap.add_argument("--tx", nargs="+", required=True)
    ap.add_argument("--k-pool", type=int, default=100)
    ap.add_argument("--hull-labels", nargs="+")
    ap.add_argument("--min-cluster-cells", type=int, default=30)
    ap.add_argument("--min-gene-umi", type=int, default=20)
    ap.add_argument("--n-jobs", type=int, default=8)
    o = ap.parse_args()

    G = ad.read_h5ad(o.embedding)
    xy = G.obsm["X_umap"]; lab = G.obs["label"].astype(str).values
    pos = {b: i for i, b in enumerate(G.obs_names)}
    emb_donors = set(G.obs["donor"].astype(str))
    n = G.n_obs
    depth = np.zeros(n); on = np.zeros(n); off = np.zeros(n)

    PD, probes_of, ex_donors = load_panel_a(o)
    panel_b_probes, seen = None, 0
    for d, path in probe_paths(o.probe_root, o.dataset, o.samples):
        a = ad.read_h5ad(path)
        lb = read_labels(o.label_root, o.dataset, d, o.drop_doublets)
        keep = [b for b in a.obs_names if b in lb.index]
        a = a[keep]
        labs = lb.loc[keep].values
        X = a.X.tocsr()
        gn = a.var["gene_name"].astype(str).values
        names = np.array(a.var_names)
        # panel b: embedding donors only
        if d in emb_donors:
            tx = a.var["tx_set"].astype(object).fillna("").astype(str).values
            gsel = gn == o.gene
            reads = gsel & np.array([all(t in s.split(",") for t in o.tx) for s in tx])
            if panel_b_probes is None:
                panel_b_probes = {"on": [s.split("|")[-1] for s in names[reads]],
                                  "off": [s.split("|")[-1] for s in names[gsel & ~reads]]}
            idx = np.array([pos[f"{d}_{b}"] for b in keep if f"{d}_{b}" in pos])
            inemb = np.array([f"{d}_{b}" in pos for b in keep])
            depth[idx] = np.asarray(X[inemb].sum(axis=1)).ravel()
            on[idx] = np.asarray(X[inemb][:, reads].sum(axis=1)).ravel()
            off[idx] = np.asarray(X[inemb][:, gsel & ~reads].sum(axis=1)).ravel()
            seen += len(idx)
        print(f"{d}: {len(keep)} labelled cells", flush=True)
    if seen != n:
        raise SystemExit(f"embedding has {n} cells, the probe matrices supplied {seen}")
    if not panel_b_probes["on"]:
        raise SystemExit(f"no {o.gene} probe reads all of {o.tx}")

    nn = NearestNeighbors(n_neighbors=o.k_pool, n_jobs=o.n_jobs).fit(G.obsm["X_pca"])
    _, ind = nn.kneighbors(G.obsm["X_pca"])
    C = sp.csr_matrix((np.ones(ind.size), (np.repeat(np.arange(n), o.k_pool), ind.ravel())),
                      shape=(n, n))
    pool = lambda v: np.asarray(C @ v.reshape(-1, 1)).ravel()
    pdepth = pool(depth)
    P4on = pool(on) / pdepth * 1e4
    P4off = pool(off) / pdepth * 1e4
    m = lab == o.label
    check_donors(ex_donors, o.dataset, emb_donors)
    panel_a = {}
    for (gene, pid), v in sorted(PD.items()):
        ins = np.array([t[0] for t in v]); out = np.array([t[1] for t in v])
        panel_a[f"{gene}|{pid}"] = {"donors": len(v), "share_label": float(ins.mean()),
                                    "share_rest": float(out.mean()),
                                    "donors_higher_in_label": int((ins > out).sum())}
    summary = {"n_cells": int(n), "donors": sorted(emb_donors), "k_pool": o.k_pool,
               "panel_a_examples": o.panel_a,
               "panel_a_donors": sorted(ex_donors[o.panel_a[0]]), "panel_a_donors_by_example": ex_donors,
               "panel_a": panel_a,
               "panel_b_probes": panel_b_probes,
               "panel_b_on_cp10k_label": float(P4on[m].mean()),
               "panel_b_on_cp10k_rest": float(P4on[~m].mean()),
               "panel_b_off_cp10k_label": float(P4off[m].mean()),
               "panel_b_off_cp10k_rest": float(P4off[~m].mean()),
               "label_counts": {k: int(v) for k, v in
                                zip(*np.unique(lab, return_counts=True))}}
    for (gene, pid), v in sorted(PD.items()):
        ins = np.array([t[0] for t in v]); out = np.array([t[1] for t in v])
        print(f"{gene:7s} {pid}  donors={len(ins)}  cell type {ins.mean():.3f} / other "
              f"{out.mean():.3f}  higher in {(ins > out).sum()}")
    print(f"panel b {o.gene} probes reading {o.tx}: {panel_b_probes['on']}; others "
          f"{panel_b_probes['off']}")
    print(f"panel b cp10k in {o.label} {summary['panel_b_on_cp10k_label']:.3f} / rest "
          f"{summary['panel_b_on_cp10k_rest']:.3f}")

    hulls = o.hull_labels or [o.label]
    H = {cl: island(xy, lab, cl) for cl in hulls}
    NOTE = {cl: cl.replace(" cells", "").replace("Naive", "naive").replace("Memory", "memory")
            for cl in hulls}
    GROUPS = [(ex.split("|")[2], ex.split("|")[1], probes_of[ex.split("|")[2]]) for ex in o.panel_a]
    # one sub-panel per example: shares per donor + the DE call per probe
    import json
    PANEL_A = []
    for ex, rj in zip(o.panel_a, o.panel_a_json):
        role, path = rj.split("=", 1)
        ds, label, gene = ex.split("|")
        st = json.load(open(path))["genes"][gene]
        info = {p["probe"]: p for p in st["probes"]}
        PANEL_A.append({"role": role, "dataset": ds, "label": label, "gene": gene,
                        "gene_fc": st["gene_fold_change"], "gene_padj": st["gene_padj"],
                        "probes": [{"probe": pid, "in": [t[0] for t in PD[(gene, pid)]],
                                    "out": [t[1] for t in PD[(gene, pid)]],
                                    "fc": info[pid]["fold_change"] if info[pid]["fold_change"] is not None else float("nan"),
                                    "sig": bool(info[pid]["sig"]),
                                    "floor": bool(info[pid]["passes_floor"]),
                                    "tx": f"{info[pid]['n_tx_hit']}/{info[pid]['n_tx_gene']}"}
                                   for pid in probes_of[gene] if (gene, pid) in PD]})
    with open(o.out, "wb") as fh:
        pickle.dump(dict(xy=xy, lab=lab, P4on=P4on, P4off=P4off, PD=PD, GROUPS=GROUPS, PANEL_A=PANEL_A,
                         H=H, NOTE=NOTE, label=o.label, gene=o.gene, summary=summary),
                    fh, protocol=4)
    import json
    json.dump(summary, open(os.path.splitext(o.out)[0] + ".summary.json", "w"), indent=1)
    print(f"wrote {o.out}", flush=True)


if __name__ == "__main__":
    main()
