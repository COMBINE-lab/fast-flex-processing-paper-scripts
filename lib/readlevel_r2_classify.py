#!/usr/bin/env python3
"""Step 13 (c): cross-tabulate Cell Ranger's idealized first-stage probe rule against the real piscem
placement of the same reads.

usage: readlevel_r2_classify.py <out.json> <prefix>...
  each <prefix> is one sample: P.cNN.{json,univ.npz} from lib/mechanism_readlevel.py --chunk NN (all chunks
  of the BAM), P_map/ (one piscem map-sc over every chunk's synthetic pairs: map_info.json,
  unmapped_bc_count.bin) and P_pairs.{bin,genes.txt} (rad_dump map on P_map/map.rad).  The synthetic R1
  cell barcode of a read is (chunk << local_bits) | its ordinal within the chunk.

Per universe read:
  cr_keep    both 25-nt halves within Hamming 1 of CR's assigned probe (ungapped, idealized stage 1)
  sa_keep    piscem placed the read and CR's assigned gene is among the genes it was placed on
  sa status  placed (same gene / other gene only), unplaced (in unmapped_bc_count.bin), or
             poisoned by a decoy (in neither the RAD nor the unmapped ledger)
  proxy_keep the exact-k-mer rule (comparison only)
"""
import json
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from readlevel_ledger import iter_unmapped  # noqa: E402

PAIR = np.dtype([("c", "<u8"), ("g", "<u4")])


def load_chunks(prefix):
    first = json.load(open(f"{prefix}.c00.json"))
    metas = [json.load(open(f"{prefix}.c{i:02d}.json")) for i in range(first["nchunks"])]
    for i, m in enumerate(metas):
        if (m["chunk"], m["nchunks"], m["local_bits"], m["genes"], m["k"]) != \
           (i, first["nchunks"], first["local_bits"], first["genes"], first["k"]):
            raise SystemExit(f"{prefix}.c{i:02d}: chunk metadata disagrees with chunk 00")
    counters = {}
    for m in metas:
        for k, v in m["counters"].items():
            counters[k] = counters.get(k, 0) + v
    return first, metas, counters


def classify(prefix):
    meta, metas, counters = load_chunks(prefix)
    zs = [np.load(f"{prefix}.c{i:02d}.univ.npz") for i in range(len(metas))]
    sizes = np.array([z["lhs_sub"].size for z in zs], dtype=np.int64)
    if any(int(sz) != m["counters"].get("universe", 0) for sz, m in zip(sizes, metas)):
        raise SystemExit(f"{prefix}: a chunk's universe rows != its counter")
    offs = np.concatenate([[0], np.cumsum(sizes)[:-1]])
    lb = meta["local_bits"]
    U = int(sizes.sum())
    lhs, rhs, gene = (np.concatenate([z[k] for z in zs]) for k in ("lhs_sub", "rhs_sub", "gene"))
    prox = np.concatenate([z["proxy_keep"] for z in zs]).astype(bool)
    del zs

    def to_global(code):
        ch, loc = code >> lb, code & ((1 << lb) - 1)
        if ch.size and (ch.max() >= sizes.size or (loc >= sizes[ch]).any()):
            raise SystemExit(f"{prefix}: a synthetic barcode decodes outside its chunk")
        return offs[ch] + loc
    mi = json.load(open(f"{prefix}_map/map_info.json"))
    if mi["num_reads"] != U:
        raise SystemExit(f"{prefix}: piscem saw {mi['num_reads']} reads, universe is {U}")

    ugenes = {g: i for i, g in enumerate(meta["genes"])}
    rgenes = [l.strip() for l in open(f"{prefix}_pairs.genes.txt") if l.strip()]
    r2u = np.array([ugenes[g] for g in rgenes], dtype=np.int64)   # every RAD gene is a reference gene
    pairs = np.fromfile(f"{prefix}_pairs.bin", dtype=PAIR)
    c = to_global(pairs["c"].astype(np.int64))
    placed = np.zeros(U, bool)
    placed[c] = True
    same = np.zeros(U, bool)
    hit = r2u[pairs["g"]] == gene[c]
    same[c[hit]] = True
    del pairs, c, hit

    unplaced = np.zeros(U, bool)
    n_led = 0
    for bcs, cnt in iter_unmapped(f"{prefix}_map/unmapped_bc_count.bin"):
        o = int(to_global(np.array([bcs[-1]], dtype=np.int64))[0])
        if cnt != 1 or unplaced[o]:
            raise SystemExit(f"{prefix}: barcode {bcs[-1]} appears more than once in the unmapped ledger")
        unplaced[o] = True
        n_led += 1
    if int(placed.sum()) != mi["num_mapped"]:
        raise SystemExit(f"{prefix}: {int(placed.sum())} placed ordinals != num_mapped {mi['num_mapped']}")
    if (placed & unplaced).any():
        raise SystemExit(f"{prefix}: an ordinal is both placed and unplaced")
    poisoned = ~(placed | unplaced)
    if int(poisoned.sum()) != mi["num_poisoned"]:
        raise SystemExit(f"{prefix}: {int(poisoned.sum())} absent ordinals != num_poisoned {mi['num_poisoned']} "
                         "(a synthetic R1 failed to parse)")

    cr = (lhs <= 1) & (rhs <= 1)
    both_gt1 = (lhs > 1) & (rhs > 1)
    one_exact_other_gt1 = (np.minimum(lhs, rhs) == 0) & (np.maximum(lhs, rhs) > 1)
    s_only = ~cr & same
    c_only = cr & ~same
    p_s_only = ~cr & prox
    p_c_only = cr & ~prox
    n = lambda m: int(m.sum())  # noqa: E731
    return {
        "universe": U,
        "bam_counters": counters, "bam_chunks": len(metas),
        "piscem_num_mapped": mi["num_mapped"], "piscem_num_poisoned": mi["num_poisoned"],
        "piscem_unmapped_ledger": n_led,
        "agree_keep": n(cr & same), "agree_reject": n(~cr & ~same),
        "simpleaf_only": n(s_only),
        "simpleaf_only_both_halves_gt1_sub": n(s_only & both_gt1),
        "simpleaf_only_one_half_exact_other_gt1_sub": n(s_only & one_exact_other_gt1),
        "simpleaf_only_other": n(s_only & ~both_gt1 & ~one_exact_other_gt1),
        "cr_only": n(c_only),
        "cr_only_unplaced": n(c_only & unplaced),
        "cr_only_poisoned": n(c_only & poisoned),
        "cr_only_other_gene_only": n(c_only & placed),
        "cr_only_halves_1_1_sub": n(c_only & (lhs == 1) & (rhs == 1)),
        "cr_keep": n(cr), "simpleaf_keep": n(same),
        "simpleaf_placed_any_gene": n(placed), "simpleaf_unplaced": n(unplaced), "simpleaf_poisoned": n(poisoned),
        "proxy": {
            "simpleaf_only": n(p_s_only), "cr_only": n(p_c_only),
            "proxy_keep_but_piscem_not_same_gene": n(prox & ~same),
            "piscem_same_gene_but_proxy_reject": n(same & ~prox),
        },
    }


def main():
    out, prefixes = sys.argv[1], sys.argv[2:]
    per = {p.rsplit("/", 1)[-1]: classify(p) for p in prefixes}
    tot = {}
    for d in per.values():
        for k, v in d.items():
            if isinstance(v, int):
                tot[k] = tot.get(k, 0) + v
            elif k == "proxy":
                tp = tot.setdefault("proxy", {})
                for kk, vv in v.items():
                    tp[kk] = tp.get(kk, 0) + vv
    if tot["agree_keep"] + tot["agree_reject"] + tot["simpleaf_only"] + tot["cr_only"] != tot["universe"]:
        raise SystemExit("R2 cross-tab does not partition the universe")
    if tot["cr_only_unplaced"] + tot["cr_only_poisoned"] + tot["cr_only_other_gene_only"] != tot["cr_only"]:
        raise SystemExit("cr_only sub-classes do not sum")
    ks = {json.load(open(f"{p}.c00.json"))["k"] for p in prefixes}
    if len(ks) != 1:
        raise SystemExit(f"samples were scanned with different k: {ks}")
    json.dump({"k": ks.pop(), "total": tot, "per_sample": per}, open(out, "w"), indent=2)
    print(json.dumps(tot, indent=2))


if __name__ == "__main__":
    main()
