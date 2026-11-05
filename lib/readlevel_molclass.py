#!/usr/bin/env python3
"""Step 13 (d): molecule-level gene reassignment, simpleaf vs Cell Ranger totals.

usage: readlevel_molclass.py <af_mol prefix> <sample_info.json> <out.json> <sample>=<bam pass prefix>,<barcodes.csv> ...
  (<bam pass prefix>.cNN.{json,crmol.npz} for every chunk NN of that sample's BAM)

A molecule is a (sample, cell, raw UMI) key, restricted to Cell Ranger's called cells of that sample.
  simpleaf  rad_dump mol on the final collated RAD: one gene per key by the alevin-fry cr-like vote,
            or AMBIG on a tie.  b0 in the collated RAD is the sample's index in sample_info.json.
  CR        every distinct (key, gene) of the reads Cell Ranger confidently assigned (fx tag), from
            lib/mechanism_readlevel.py; reads whose gene is absent from the simpleaf reference are
            counted separately; they cannot be compared.
CR-side classes of each distinct CR (key, gene):  recovered (simpleaf key present, same gene),
reassigned_out (present, a different gene), lost_tie (present, AMBIG), lost_absent (key absent).
simpleaf-side classes of each non-AMBIG simpleaf molecule:  recovered, reassigned_in (key is a CR
molecule but CR never gave it this gene), novel (key not a CR molecule).
"""
import json
import sys

import numpy as np

REC = np.dtype([("k", "<u8"), ("g", "<u4")])
AMBIG = np.uint32(0xFFFFFFFF)
MASK56 = np.uint64((1 << 56) - 1)
ENC = {"A": 0, "C": 1, "G": 2, "T": 3}


def code_of(s):
    c = 0
    for ch in s:
        c = (c << 2) | ENC[ch]
    return c


def cells_of(csv):
    out = set()
    for line in open(csv):
        bc = line.strip().split(",")[-1]
        if bc and bc != "barcode":
            out.add(code_of(bc.split("-")[0][:16]))
    return np.array(sorted(out), dtype=np.uint64)


def main():
    af_prefix, sinfo, out = sys.argv[1:4]
    specs = dict(a.split("=", 1) for a in sys.argv[4:])
    si = json.load(open(sinfo))
    names = [s["name"] for s in si["samples"]]
    bmap = {}
    with open(f"{af_prefix}.bmap.tsv") as f:
        next(f)
        for line in f:
            b0, bidx, nrec = (int(x) for x in line.split("\t"))
            if b0 >= len(names) or si["samples"][b0]["collatable_reads"] != nrec:
                raise SystemExit(f"b0 {b0}: {nrec} records do not match sample_info collatable reads")
            bmap[bidx] = b0
    genes = [l.strip() for l in open(f"{af_prefix}.genes.txt") if l.strip()]
    gidx = {g: i for i, g in enumerate(genes)}
    G = np.uint64(len(genes))
    af = np.fromfile(f"{af_prefix}.bin", dtype=REC)
    af_sample = np.array([bmap[b] for b in range(len(bmap))], dtype=np.int64)[(af["k"] >> np.uint64(56)).astype(np.int64)]

    per = {}
    for name, spec in specs.items():
        prefix, bcsv = spec.split(",")
        si_idx = names.index(name)
        cells = cells_of(bcsv)
        sel = af_sample == si_idx
        ak, ag = af["k"][sel] & MASK56, af["g"][sel]
        if np.unique(ak).size != ak.size:
            raise SystemExit(f"{name}: simpleaf molecule keys are not unique")
        in_cells = np.isin(ak >> np.uint64(24), cells)
        off_axis = int((~in_cells).sum())
        ak, ag = ak[in_cells], ag[in_cells]

        meta = json.load(open(f"{prefix}.c00.json"))
        zs = [np.load(f"{prefix}.c{i:02d}.crmol.npz") for i in range(meta["nchunks"])]
        for i in range(meta["nchunks"]):
            if json.load(open(f"{prefix}.c{i:02d}.json"))["features"] != meta["features"]:
                raise SystemExit(f"{prefix}.c{i:02d}: features list differs from chunk 00")
        feat2g = np.array([gidx.get(f, -1) for f in meta["features"]], dtype=np.int64)
        cg = feat2g[np.concatenate([z["feature"] for z in zs])]
        absent = cg < 0
        ck = np.concatenate([z["key"] for z in zs])[~absent]
        cg = cg[~absent].astype(np.uint64)
        del zs
        if not np.isin(ck >> np.uint64(24), cells).all():
            raise SystemExit(f"{name}: a CR molecule lies outside CR's called cells")

        all_k = np.unique(np.concatenate([ak, ck]))
        a_id = np.searchsorted(all_k, ak).astype(np.uint64)
        c_comp = np.unique(np.searchsorted(all_k, ck).astype(np.uint64) * G + cg)
        cu_id, cu_g = c_comp // G, (c_comp % G).astype(np.uint32)
        af_by_id = np.full(all_k.size, 0xFFFFFFFE, dtype=np.uint32)   # 0xFFFFFFFE = key absent in simpleaf
        af_by_id[a_id] = ag
        here = af_by_id[cu_id]
        m_absent = here == np.uint32(0xFFFFFFFE)
        m_tie = here == AMBIG
        m_same = here == cu_g
        m_reas = ~m_absent & ~m_tie & ~m_same

        res = ag != AMBIG
        a_comp = a_id[res] * G + ag[res].astype(np.uint64)
        in_cr = np.isin(a_comp, c_comp)
        key_in_cr = np.isin(a_id[res][~in_cr], np.unique(cu_id))
        row = {
            "cr_called_cells": int(cells.size),
            "cr_reads_gene_absent_in_simpleaf_reference": int(absent.sum()),
            "cr_molecules": int(c_comp.size),
            "cr_multigene_surplus": int(c_comp.size - np.unique(cu_id).size),
            "simpleaf_molecules": int(res.sum()),
            "simpleaf_ambiguous_molecules": int((~res).sum()),
            "simpleaf_molecules_off_axis_dropped": off_axis,
            "recovered": int(m_same.sum()),
            "reassigned_out": int(m_reas.sum()),
            "lost_tie": int(m_tie.sum()),
            "lost_absent": int(m_absent.sum()),
            "reassigned_in": int(key_in_cr.sum()),
            "novel": int((~key_in_cr).sum()),
        }
        if int(in_cr.sum()) != row["recovered"]:
            raise SystemExit(f"{name}: recovered differs between sides ({int(in_cr.sum())} vs {row['recovered']})")
        if row["recovered"] + row["reassigned_out"] + row["lost_tie"] + row["lost_absent"] != row["cr_molecules"]:
            raise SystemExit(f"{name}: CR partition open")
        if row["recovered"] + row["reassigned_in"] + row["novel"] != row["simpleaf_molecules"]:
            raise SystemExit(f"{name}: simpleaf partition open")
        per[name] = row
        print(name, json.dumps(row), flush=True)

    tot = {k: sum(r[k] for r in per.values()) for k in next(iter(per.values()))}
    json.dump({"total": tot, "per_sample": per}, open(out, "w"), indent=2)
    print(json.dumps(tot, indent=2))


if __name__ == "__main__":
    main()
