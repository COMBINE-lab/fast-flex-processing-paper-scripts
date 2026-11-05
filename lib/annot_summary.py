#!/usr/bin/env python3
"""Methods "Probe annotation" numbers, and the gate that the probe matrices carry this annotation.

From the vendor probe-set CSV and the relocation table written by probe_annotate.py:
  region counts (`spliced` / `unspliced`) over all probes and over `included=TRUE` probes
  (Methods quotes all 54,580 human v1.1.0 probes: 7,605 / 46,975), genomic_class counts, the
  region x genomic_class cross-tab (the two groupings are not interchangeable), and the number of
  isoform-informative genes (probes spanning >= 2 distinct isoform sets).
Gate (--h5ad): every probe of the matrix must carry the same genomic_class and tx_set as the
table; a mismatch exits non-zero, because the analyses would then describe another annotation.

Usage: annot_summary.py --probe-csv CSV --annot TSV --out JSON [--h5ad H5AD]
"""
import argparse
import json

import pandas as pd


def read_probe_csv(path):
    return pd.read_csv(path, comment="#")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe-csv", required=True); ap.add_argument("--annot", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--h5ad")
    o = ap.parse_args()

    csv = read_probe_csv(o.probe_csv)
    ann = pd.read_csv(o.annot, sep="\t").set_index("probe_id")
    ann["tx_set"] = ann["tx_set"].fillna("")
    inc = csv["included"].astype(str).str.upper() == "TRUE"
    a_inc = ann[ann.included.astype(str).str.upper() == "TRUE"]
    ntx = a_inc.groupby("gene_id").tx_set.nunique()
    res = {
        "probe_csv": o.probe_csv, "n_probes": int(len(csv)), "n_included": int(inc.sum()),
        "region_all": {k: int(v) for k, v in csv["region"].value_counts().items()},
        "region_included": {k: int(v) for k, v in csv.loc[inc, "region"].value_counts().items()},
        "genomic_class_all": {k: int(v) for k, v in ann.genomic_class.value_counts().items()},
        "genomic_class_included": {k: int(v) for k, v in
                                   a_inc.genomic_class.value_counts().items()},
        "region_x_class_included": {f"{r}|{c}": int(v) for (r, c), v in
                                    a_inc.groupby(["region_csv", "genomic_class"]).size().items()},
        "n_genes_included": int(a_inc.gene_id.nunique()),
        "n_isoform_informative_genes": int((ntx >= 2).sum()),
    }
    # probe_annotate.py writes only probes whose gene is in the GTF; the rest are counted here
    extra = ann.index.difference(csv["probe_id"])
    if len(extra):
        raise SystemExit(f"{len(extra)} annotated probes are not in the CSV")
    res["n_annotated"] = int(len(ann))
    res["n_gene_absent_from_gtf"] = int(len(csv) - len(ann))
    res["n_included_gene_absent_from_gtf"] = int(inc.sum() - len(a_inc))
    if o.h5ad:
        import anndata as ad
        v = ad.read_h5ad(o.h5ad, backed="r").var.copy()
        v["tx_set"] = v["tx_set"].astype(object).fillna("").astype(str)
        miss = v.index.difference(ann.index)
        if len(miss):
            raise SystemExit(f"{len(miss)} matrix probes absent from the annotation, e.g. "
                             f"{list(miss[:3])}")
        a = ann.loc[v.index]
        bad = (a.genomic_class.astype(str).values != v.genomic_class.astype(str).values) | \
              (a.tx_set.astype(str).values != v.tx_set.values)
        if bad.any():
            raise SystemExit(f"{int(bad.sum())} of {len(v)} matrix probes carry another "
                             f"genomic_class/tx_set than {o.annot}")
        res["h5ad_gate"] = {"h5ad": o.h5ad, "n_probes": int(len(v)), "mismatches": 0}
    json.dump(res, open(o.out, "w"), indent=1)
    print(json.dumps({k: res[k] for k in ("n_probes", "n_included", "region_all",
                                          "region_included", "genomic_class_included")}))


if __name__ == "__main__":
    main()
