#!/usr/bin/env python3
"""
CellTypist annotation of the Flex gene matrices, harmonised across samples.

The August campaign clustered every (dataset, sample) with `scbp_pipeline.py` and kept the
per-cell obs table (QC-passing cells, scran size factors, Leiden at three resolutions).
Leiden labels are per-sample and therefore not comparable between samples, which is a problem
because the probe-level tests need cell groups that are the same thing in every sample and use
samples as replicates. CellTypist labels supply exactly that, and running it in
`majority_voting` mode over the campaign's own `leiden_res1_0` keeps the label consistent
within each cluster the clustering already found.

Annotation runs on the **gene** matrix from the simpleaf quant -- the same quant the probe
matrix is derived from, so gene and probe counts are one measurement.

No CellTypist model covers mouse spleen (the mouse models are brain, gut, liver and
dendritic-cell), so mouse symbols are mapped to human orthologs via the MGI homology report
and the human immune model is used. That is CellTypist's own recommendation for cross-species
annotation, and it is recorded per run so the writeup can state it.

Usage:
  celltype_annotate.py --dataset 16k_mouse_Spleno --outdir DIR [--sample S] [--model M]
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flexio

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = flexio._bench()
CLUSTERING = os.path.join(BENCH, "concordance", "clustering")
MGI = os.path.join(os.path.dirname(HERE), "data", "HOM_MouseHumanSequence.rpt")   # = config.env MGI_HOMOLOGY

# Immune_All_Low is the only model spanning the immune compartments these datasets contain.
# 40k_DTC also holds tumour and stromal cells, which it cannot represent -- those come back
# with low confidence and must be resolved by marker inspection, not by trusting the label.
DEFAULT_MODEL = "Immune_All_Low.pkl"
ORGANISM = {"16k_mouse_Spleno": "mouse"}


def obs_table(dataset, sample, tool="simpleaf"):
    """The campaign's post-QC obs table for one sample."""
    stem = f"{dataset}__{sample}__{tool}__common"
    tsv = os.path.join(CLUSTERING, stem + ".obs.tsv.gz")
    if os.path.exists(tsv):
        return pd.read_csv(tsv, sep="\t", index_col=0)
    h5 = os.path.join(CLUSTERING, stem + ".h5ad")
    if os.path.exists(h5):
        import anndata as ad
        return ad.read_h5ad(h5, backed="r").obs.copy()
    return None


def mouse_to_human():
    """MGI homology report -> {mouse symbol: human symbol} for 1:1 classes."""
    df = pd.read_csv(MGI, sep="\t", usecols=["DB Class Key", "Common Organism Name", "Symbol"])
    df["org"] = np.where(df["Common Organism Name"].str.contains("mouse", case=False),
                         "mouse", "human")
    m = {}
    for key, g in df.groupby("DB Class Key"):
        mm = g.loc[g.org == "mouse", "Symbol"].tolist()
        hh = g.loc[g.org == "human", "Symbol"].tolist()
        if len(mm) == 1 and len(hh) == 1:                # 1:1 orthologs only
            m[mm[0]] = hh[0]
    return m


def build_gene_adata(dataset, sample, obs, organism):
    """Post-QC cells x genes, gene symbols on the var axis, CellTypist-normalised."""
    import anndata as ad
    import scanpy as sc

    qdir = flexio.simpleaf_quant_dir(dataset)
    if not qdir:
        raise SystemExit(f"no simpleaf quant for {dataset}")
    loaded = flexio.load_alevin(qdir, sample)
    if loaded is None:
        raise SystemExit(f"no simpleaf rows for {dataset}/{sample}")
    M, bcs, gene_ids = loaded

    row = flexio.dataset_row(dataset)
    name2id = flexio.probe_name_to_id(row["probe_csv"])
    id2name = {}
    for nm, gid in name2id.items():
        id2name.setdefault(gid, nm)

    bi = {b: i for i, b in enumerate(bcs)}
    cells = [b for b in obs.index if b in bi]
    A = ad.AnnData(M[[bi[b] for b in cells], :].tocsr().astype(np.float32))
    A.obs_names = cells
    A.var_names = gene_ids
    A.obs = A.obs.join(obs.loc[cells])
    A.var["symbol"] = [id2name.get(g, g) for g in gene_ids]

    if organism == "mouse":
        conv = mouse_to_human()
        A.var["symbol"] = [conv.get(s, "") for s in A.var["symbol"]]
        A = A[:, A.var["symbol"] != ""].copy()

    # collapse to symbols (CellTypist keys on symbols, and duplicates must be summed)
    A.var_names = A.var["symbol"].values
    A.var_names_make_unique()
    A.layers["counts"] = A.X.copy()
    sc.pp.normalize_total(A, target_sum=1e4)             # CellTypist's required input scale
    sc.pp.log1p(A)
    return A


def annotate_sample(dataset, sample, model, outdir):
    import celltypist

    organism = ORGANISM.get(dataset, "human")
    obs = obs_table(dataset, sample)
    if obs is None:
        print(f"[annot] {dataset}/{sample}: no campaign obs table, skipped", flush=True)
        return None
    A = build_gene_adata(dataset, sample, obs, organism)

    over = "leiden_res1_0" if "leiden_res1_0" in A.obs.columns else None
    if over is not None:
        A.obs[over] = A.obs[over].astype(str)
    pred = celltypist.annotate(A, model=model, majority_voting=True, over_clustering=over)
    out = pred.predicted_labels.copy()
    out.index.name = "barcode"
    out["sample"] = sample
    for c in ("leiden_res0_25", "leiden_res0_5", "leiden_res1_0", "size_factors",
              "total_counts", "pct_counts_mt", "scDblFinder_class"):
        if c in A.obs.columns:
            out[c] = A.obs[c].values

    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{sample}.celltypist.tsv.gz")
    out.to_csv(path, sep="\t")
    col = "majority_voting" if "majority_voting" in out.columns else "predicted_labels"
    print(f"[annot] {dataset}/{sample}: {len(out)} cells, {out[col].nunique()} labels, "
          f"organism={organism}, genes_used={A.n_vars}", flush=True)
    return {"path": path, "n_cells": int(len(out)), "n_labels": int(out[col].nunique()),
            "organism": organism, "model": model, "genes_used": int(A.n_vars),
            "over_clustering": over,
            "top_labels": out[col].value_counts().head(15).to_dict()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--sample")
    ap.add_argument("--model", default=DEFAULT_MODEL)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    samples = [a.sample] if a.sample else [s for s, _ in flexio.samples_of(a.dataset)]
    res = {}
    for s in samples:
        r = annotate_sample(a.dataset, s, a.model, a.outdir)
        if r:
            res[s] = r
    with open(os.path.join(a.outdir, "annotate.json"), "w") as f:
        json.dump(res, f, indent=2)
    sys.exit(0 if res else 1)


if __name__ == "__main__":
    main()
