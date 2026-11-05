#!/usr/bin/env python3
"""Download the CellTypist model into $CELLTYPIST_FOLDER and record what the Methods "Cell-type annotation" state about it.

Writes JSON: celltypist version, model file + sha256 + models.json entry (version/date), number of model
features and cell types, and the mouse cross-species numbers: MGI homology classes with both species, 1:1
classes, mouse probe-panel genes with a 1:1 human ortholog, and how many of those are model features.

Usage: celltypist_model_check.py --model Immune_All_Low.pkl --mgi RPT --mouse-probe-csv CSV --out JSON
"""
import argparse, hashlib, json, os, sys

import numpy as np
import pandas as pd
import celltypist
from celltypist import models


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True); ap.add_argument("--mgi", required=True)
    ap.add_argument("--mouse-probe-csv", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if "CELLTYPIST_FOLDER" not in os.environ:
        raise SystemExit("set CELLTYPIST_FOLDER (the durable model home)")
    path = os.path.join(models.models_path, a.model)
    if not os.path.exists(path):
        models.download_models(model=a.model, force_update=True)
    sha = hashlib.sha256(open(path, "rb").read()).hexdigest()
    meta = json.load(open(os.path.join(models.models_path, "models.json")))
    entry = next((m for m in meta.get("models", []) if m.get("filename") == a.model), None)
    m = models.Model.load(a.model)
    feats = set(m.features)

    df = pd.read_csv(a.mgi, sep="\t", usecols=["DB Class Key", "Common Organism Name", "Symbol"])
    df["org"] = np.where(df["Common Organism Name"].str.contains("mouse", case=False), "mouse", "human")
    both = one = 0
    m2h = {}
    for _, g in df.groupby("DB Class Key"):
        mm, hh = g.loc[g.org == "mouse", "Symbol"].tolist(), g.loc[g.org == "human", "Symbol"].tolist()
        if mm and hh:
            both += 1
            if len(mm) == 1 and len(hh) == 1:
                one += 1
                m2h[mm[0]] = hh[0]
    rows = [l for l in open(a.mouse_probe_csv) if not l.startswith("#")]
    hdr = rows[0].rstrip("\n").split(",")
    gi = hdr.index("gene_name")
    panel = {l.rstrip("\n").split(",")[gi] for l in rows[1:] if l.strip()}
    mapped = {m2h[g] for g in panel if g in m2h}
    out = {"celltypist_version": celltypist.__version__, "model": a.model, "model_path": path, "sha256": sha,
           "models_json_entry": entry, "n_features": len(feats), "n_cell_types": len(m.cell_types),
           "mgi_classes_with_both_species": both, "mgi_one_to_one_classes": one,
           "mouse_panel_genes": len(panel), "mouse_panel_genes_with_1to1_ortholog": len([g for g in panel if g in m2h]),
           "mapped_orthologs_in_model_features": len(mapped & feats)}
    out["mouse_panel_1to1_pct"] = 100.0 * out["mouse_panel_genes_with_1to1_ortholog"] / len(panel)
    out["model_features_reached_pct"] = 100.0 * out["mapped_orthologs_in_model_features"] / len(feats)
    json.dump(out, open(a.out, "w"), indent=2)
    print(json.dumps({k: v for k, v in out.items() if k != "models_json_entry"}, indent=2))


if __name__ == "__main__":
    main()
