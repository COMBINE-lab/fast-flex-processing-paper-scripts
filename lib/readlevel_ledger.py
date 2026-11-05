#!/usr/bin/env python3
"""Step 13 helpers: geometry parsing, unmapped-barcode ledger, and the R1 stage-wise read ledger.

  readlevel_ledger.py geom '<GEOM_V2>'
      print shell assignments CB_LEN UMI_LEN SP_LO SP_HI ANCHOR MAXD SBC_LEN (fails on any other shape)
  readlevel_ledger.py ledger <map_info.json> <unmapped_bc_count.bin> <af_quant dir> <out.json>
      read-level ledger of the final run (map_info.json + unmapped_bc_count.bin + sample_info.json)
  readlevel_ledger.py reconcile <ledger.json> <scan part.json>... <out.json>
      merge the R1 scanner parts and require their drop to equal the ledger's R1 drop exactly
"""
import json
import re
import struct
import sys

GEOM_RE = re.compile(
    r"^1\{b\[(\d+)\]u\[(\d+)\]x\[(\d+)-(\d+)\]hamming\(f\[([ACGT]+)\],(\d+)\)s\[(\d+)\]x:\}2\{r:\}$")
TYPE_WIDTH = {1: 1, 2: 2, 3: 4, 4: 8}   # piscem-rs unmapped_bc_count.bin type ids -> bytes


def parse_geom(geom):
    m = GEOM_RE.match(geom)
    if not m:
        raise SystemExit(f"GEOM_V2 is not of the form 1{{b[]u[]x[a-b]hamming(f[],d)s[]x:}}2{{r:}}: {geom}")
    cb, umi, lo, hi, anchor, maxd, sbc = m.groups()
    return dict(CB_LEN=int(cb), UMI_LEN=int(umi), SP_LO=int(lo), SP_HI=int(hi),
                ANCHOR=anchor, MAXD=int(maxd), SBC_LEN=int(sbc))


def read_unmapped(path):
    """Sum of counts and number of records in a piscem-rs v1 unmapped_bc_count.bin."""
    with open(path, "rb") as f:
        data = f.read()
    if data[0] != 1:
        raise SystemExit(f"{path}: unmapped_bc_count version {data[0]} != 1")
    nf = data[1]
    widths = [TYPE_WIDTH[t] for t in data[2:2 + nf]]
    rec = sum(widths) + 4
    body = memoryview(data)[2 + nf:]
    if len(body) % rec:
        raise SystemExit(f"{path}: body {len(body)} bytes is not a multiple of record size {rec}")
    total = 0
    off = sum(widths)
    for i in range(0, len(body), rec):
        total += struct.unpack_from("<I", body, i + off)[0]
    return total, len(body) // rec, widths


def iter_unmapped(path):
    """Yield (barcode tuple, count) records of a piscem-rs v1 unmapped_bc_count.bin."""
    with open(path, "rb") as f:
        data = f.read()
    nf = data[1]
    widths = [TYPE_WIDTH[t] for t in data[2:2 + nf]]
    fmt = "<" + "".join({1: "B", 2: "H", 4: "I", 8: "Q"}[w] for w in widths) + "I"
    rec = struct.calcsize(fmt)
    for vals in struct.iter_unpack(fmt, memoryview(data)[2 + nf:]):
        yield vals[:-1], vals[-1]
    assert rec == sum(widths) + 4


def ledger(map_info, unmapped, af_quant, out):
    mi = json.load(open(map_info))
    si = json.load(open(f"{af_quant}/sample_info.json"))
    led, nrec, widths = read_unmapped(unmapped)
    n, mapped, poisoned = mi["num_reads"], mi["num_mapped"], mi["num_poisoned"]
    r1_drop = n - mapped - poisoned - led
    if r1_drop < 0:
        raise SystemExit(f"negative R1 drop {r1_drop}")
    sr = si["sample_routing"]
    sbc_rejected = sr["immediate_rejected"] + sr["deferred_rejected"]
    sbc_exact = sr["exact"]
    sbc_corrected = sr["structurally_unique"] + sr["deferred_accepted"]
    if sbc_exact + sbc_corrected + sbc_rejected != mapped:
        raise SystemExit(f"sample routing {sbc_exact}+{sbc_corrected}+{sbc_rejected} != num_mapped {mapped}")
    if si["matched_reads"] != sbc_exact + sbc_corrected:
        raise SystemExit(f"matched_reads {si['matched_reads']} != exact+corrected {sbc_exact + sbc_corrected}")
    per = {}
    tot = dict(cell_exact=0, cell_corrected=0, cell_not_found=0, cell_ambiguous=0, collatable=0)
    for s in si["samples"]:
        cc = s["cell_correction"]
        if cc["exact_reads"] + cc["corrected_reads"] != s["collatable_reads"]:
            raise SystemExit(f"{s['name']}: exact+corrected != collatable_reads")
        row = dict(cell_exact=cc["exact_reads"], cell_corrected=cc["corrected_reads"],
                   cell_not_found=cc["not_found_reads"], cell_ambiguous=cc["ambiguous_reads"],
                   collatable=s["collatable_reads"], num_cells=s["num_cells"])
        per[s["name"]] = row
        for k in tot:
            tot[k] += row[k]
    routed = tot["collatable"] + tot["cell_not_found"] + tot["cell_ambiguous"]
    if routed != sbc_exact + sbc_corrected:
        raise SystemExit(f"per-sample routed reads {routed} != sample-assigned reads {sbc_exact + sbc_corrected}")
    res = {
        "num_reads": n,
        "r1_parse_drop": r1_drop,
        "poisoned_by_decoy": poisoned,
        "r2_no_placement": led,
        "mapped": mapped,
        "sample_bc_rejected": sbc_rejected,
        "sample_bc_exact": sbc_exact,
        "sample_bc_corrected": sbc_corrected,
        **tot,
        "per_sample": per,
        "unmapped_bc_records": nrec,
        "unmapped_bc_field_bytes": widths,
        "mapping_rate": mapped / n,
        "r1_parse_pass_rate": (n - r1_drop) / n,
        "map_info_piscem_rs_version": mi.get("piscem_rs_version"),
        "map_info_index_sha256_seqs": mi["signatures"]["sha256_seqs"],
    }
    # the whole budget must close to zero
    closes = (n - r1_drop - poisoned - led - sbc_rejected - tot["cell_not_found"]
              - tot["cell_ambiguous"] - tot["collatable"])
    if closes != 0:
        raise SystemExit(f"read ledger does not close: residual {closes}")
    json.dump(res, open(out, "w"), indent=2)
    print(json.dumps({k: v for k, v in res.items() if k != "per_sample"}, indent=2))


def reconcile(ledger_json, parts, out):
    led = json.load(open(ledger_json))
    tot = None
    for p in parts:
        d = json.load(open(p))
        if tot is None:
            tot = d
            continue
        if d["geometry"] != tot["geometry"]:
            raise SystemExit(f"{p}: geometry differs between parts")
        tot["num_reads"] += d["num_reads"]
        for key in ("classes", "anchor_absent_subclasses", "parse_ok_by_offset_dist",
                    "anchor_absent_by_found_offset_dist"):
            for k, v in d[key].items():
                tot[key][k] = tot[key].get(k, 0) + v
        for key in ("parse_ok_sample_bc_exact_listed", "parse_ok_sample_bc_not_listed"):
            tot[key] += d[key]
        tot["parse_ok_per_sample_exact"] = [a + b for a, b in zip(tot["parse_ok_per_sample_exact"],
                                                                  d["parse_ok_per_sample_exact"])]
    if tot["num_reads"] != led["num_reads"]:
        raise SystemExit(f"scan saw {tot['num_reads']} reads, map_info has {led['num_reads']}")
    drop_scan = tot["num_reads"] - tot["classes"]["parse_ok"]
    if sum(tot["anchor_absent_subclasses"].values()) != tot["classes"]["anchor_absent_in_window"]:
        raise SystemExit("anchor-absent subclasses do not sum to the class")
    tot["reconciliation"] = {
        "formula": "r1_drop = num_reads - num_mapped - num_poisoned - sum(unmapped_bc_count.bin)",
        "r1_drop_from_map_info": led["r1_parse_drop"],
        "r1_drop_from_scan": drop_scan,
        "difference": drop_scan - led["r1_parse_drop"],
    }
    json.dump(tot, open(out, "w"), indent=2)
    print(json.dumps(tot["classes"]), json.dumps(tot["anchor_absent_subclasses"]),
          json.dumps(tot["reconciliation"]))
    if drop_scan != led["r1_parse_drop"]:
        raise SystemExit(f"R1 scan drop {drop_scan} != map_info R1 drop {led['r1_parse_drop']}")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "geom":
        for k, v in parse_geom(sys.argv[2]).items():
            print(f"{k}={v}")
    elif cmd == "ledger":
        ledger(*sys.argv[2:6])
    elif cmd == "reconcile":
        reconcile(sys.argv[2], sys.argv[3:-1], sys.argv[-1])
    else:
        raise SystemExit(f"unknown command {cmd}")
