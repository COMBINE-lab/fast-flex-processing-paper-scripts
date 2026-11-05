#!/usr/bin/env python3
"""Step 13 BAM pass over one Cell Ranger Flex per-sample BAM.

  * K is a parameter (--k).
  * The pass does not decide the simpleaf side itself.  For every universe read it writes a
    synthetic read pair (R1 = a valid Flex v2 R1 whose 16-nt cell barcode encodes the read's
    ordinal, R2 = the read in its sequenced orientation) so that the real piscem binary and index
    can place it; lib/readlevel_r2_classify.py joins piscem's placements back by ordinal.
    An exact-K-mer rule is also applied as `proxy_keep`, for comparison.
  * The same pass extracts Cell Ranger's per-read (cell16, raw UMI, gene) molecules for the gene
    reassignment analysis: an `fx` gene tag, a CB whose 16-nt core is a called cell of this
    sample, and fx in the CR features list.

Universe: primary, mapped records with a `pr` tag, a 50-nt SEQ and the probe present in
the simpleaf reference FASTA.  `o` = SEQ reverse-complemented when the record is reverse, i.e. the
read as sequenced, which is probe-sense.  CR rule (idealized first stage): each 25-nt half of `o`
within Hamming 1 of the assigned probe's half.

Outputs (prefix P):  P.json counters;  P.univ.npz (lhs_sub, rhs_sub, gene, proxy_keep per ordinal);
P.crmol.npz (key = cell16<<24 | umi, feature index);  P_R1.fastq.gz / P_R2.fastq.gz synthetic reads.
"""
import argparse
import array
import collections
import gzip
import json
import operator
import os
import subprocess
import sys
import time

import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from readlevel_ledger import parse_geom  # noqa: E402

_RC = str.maketrans("ACGTN", "TGCAN")
_ENC = {"A": 0, "C": 1, "G": 2, "T": 3}
_BASES = "ACGT"
HALF = 25


def revcomp(s):
    return s.translate(_RC)[::-1]


def hamming(a, b):
    return sum(map(operator.ne, a, b))


def kmer_ints(seq, k):
    """Yield the 2-bit-encoded k-mer ints of seq, skipping any window containing N."""
    codes = [_ENC.get(ch, -1) for ch in seq]
    mask = (1 << (2 * k)) - 1
    val, run = 0, 0
    for c in codes:
        if c < 0:
            val, run = 0, 0
            continue
        val = ((val << 2) | c) & mask
        run += 1
        if run >= k:
            yield val


def build_index(ref_fa, k):
    """probe name -> 50-nt seq; gene id -> frozenset of its probes' k-mers; ordered gene list."""
    probes = {}
    gene_kmers = collections.defaultdict(set)
    name, seq = None, []

    def flush():
        if name is None:
            return
        s = "".join(seq)
        probes[name] = s
        gene_kmers[name.split("|")[0]].update(kmer_ints(s, k))
    with open(ref_fa) as f:
        for line in f:
            if line.startswith(">"):
                flush()
                name, seq = line[1:].strip().split()[0], []
            else:
                seq.append(line.strip())
    flush()
    genes = sorted(gene_kmers)
    return probes, {g: frozenset(v) for g, v in gene_kmers.items()}, genes


def code_of(s):
    """2-bit code, first base most significant; non-ACGT -> 0."""
    c = 0
    for ch in s:
        c = (c << 2) | _ENC.get(ch, 0)
    return c


def encode_ordinal(u, n):
    return "".join(_BASES[(u >> (2 * (n - 1 - i))) & 3] for i in range(n))


def chunk_regions(bam, nchunks):
    """Deterministic split of a coordinate-sorted, indexed BAM into nchunks lists of (contig, start, end).

    Each contig is cut into equal-length pieces sized from its mapped-read count (4 pieces per chunk on
    average); consecutive pieces are packed greedily into chunks of ~equal estimated read count. A record
    belongs to the piece holding its reference_start, so the pieces partition the mapped records exactly.
    The no-coordinate records (contig None) go to the last chunk."""
    st = {s.contig: s.mapped for s in bam.get_index_statistics()}
    total = sum(st.values())
    grain = max(1, total // (4 * nchunks))
    pieces = []
    for contig, length in zip(bam.references, bam.lengths):
        m = st.get(contig, 0)
        if not m:
            continue
        k = max(1, -(-m // grain))
        step = -(-length // k)
        pieces += [(contig, a, min(a + step, length), m / k) for a in range(0, length, step)]
    chunks = [[] for _ in range(nchunks)]
    target, acc, i = total / nchunks, 0.0, 0
    for contig, a, b, est in pieces:
        if acc >= target * (i + 1) and i < nchunks - 1:
            i += 1
        chunks[i].append((contig, a, b))
        acc += est
    chunks[-1].append((None, 0, 0))
    return chunks


def iter_chunk(bam, regions):
    for contig, a, b in regions:
        if contig is None:
            yield from bam.fetch(contig="*")
            continue
        for r in bam.fetch(contig, a, b):
            if a <= r.reference_start < b:
                yield r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True, help="simpleaf probe reference ref.fa")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--geom", required=True, help="GEOM_V2, used to build the synthetic R1")
    ap.add_argument("--sbc", required=True, help="sample barcode written into every synthetic R1")
    ap.add_argument("--barcodes", required=True, help="CR sample_filtered_barcodes.csv")
    ap.add_argument("--features", required=True, help="CR features.tsv.gz (gene id in column 1)")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--limit", type=int, default=0, help="stop after N records (0 = all)")
    ap.add_argument("--chunk", type=int, default=0, help="this worker's chunk of the BAM (0-based)")
    ap.add_argument("--nchunks", type=int, default=1, help="number of coordinate chunks the BAM is split into")
    ap.add_argument("--progress", type=int, default=20_000_000)
    ap.add_argument("--gzip-threads", type=int, default=4)
    args = ap.parse_args()
    K = args.k

    g = parse_geom(args.geom)
    if len(args.sbc) != g["SBC_LEN"]:
        raise SystemExit(f"--sbc length {len(args.sbc)} != geometry s[{g['SBC_LEN']}]")
    r1_tail = "ACGT" * (g["UMI_LEN"] // 4) + "ACGT"[: g["UMI_LEN"] % 4] + "A" * g["SP_LO"] + g["ANCHOR"] + args.sbc
    r1_len = g["CB_LEN"] + len(r1_tail)
    q1 = "F" * r1_len

    t0 = time.time()
    probes, gene_kmers, genes = build_index(args.ref, K)
    gidx = {x: i for i, x in enumerate(genes)}
    print(f"[{args.bam}] k={K}: {len(probes):,} probes, {len(genes):,} genes", flush=True)

    features = [line.split("\t")[0] for line in gzip.open(args.features, "rt")]
    fidx = {x: i for i, x in enumerate(features)}
    cells = set()
    for line in open(args.barcodes):
        bc = line.strip().split(",")[-1]
        if bc and bc != "barcode":
            cells.add(bc.split("-")[0][:16])
    if not cells:
        raise SystemExit(f"no called cells in {args.barcodes}")

    gz = [subprocess.Popen(["pigz", "-1", "-p", str(args.gzip_threads), "-c"],
                           stdin=subprocess.PIPE, stdout=open(f"{args.out}_R{i}.fastq.gz", "wb"))
          for i in (1, 2)]
    w1, w2 = gz[0].stdin, gz[1].stdin
    b1, b2 = [], []

    c = collections.Counter()
    lhs_a, rhs_a, gene_a, prox_a = array.array("B"), array.array("B"), array.array("I"), array.array("B")
    mol_key, mol_feat = array.array("Q"), array.array("I")
    u = 0
    n = 0
    # synthetic cell barcode = (chunk << local_bits) | ordinal within the chunk
    chunk_bits = max(1, (args.nchunks - 1).bit_length())
    local_bits = 2 * g["CB_LEN"] - chunk_bits
    if not 0 <= args.chunk < args.nchunks or local_bits < 24:
        raise SystemExit(f"bad --chunk {args.chunk} / --nchunks {args.nchunks}")
    base = args.chunk << local_bits
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        records = bam if args.nchunks == 1 else iter_chunk(bam, chunk_regions(bam, args.nchunks)[args.chunk])
        for r in records:
            if args.limit and n >= args.limit:
                break
            n += 1
            if args.progress and n % args.progress == 0:
                print(f"  {n:,}  ({n / (time.time() - t0):,.0f} rec/s)", flush=True)
            # ---- Cell Ranger molecules (crrad_make_sam.py inclusion rule, every record) ----
            if r.has_tag("fx"):
                c["fx_reads"] += 1
                cb = r.get_tag("CB") if r.has_tag("CB") else None
                if cb and cb != "-":
                    cell = cb.split("-")[0][:16]
                    if cell in cells:
                        fi = fidx.get(r.get_tag("fx"))
                        ur = r.get_tag("UR") if r.has_tag("UR") else None
                        if fi is None:
                            c["fx_incell_not_in_features"] += 1
                        elif ur:
                            mol_key.append((code_of(cell) << 24) | code_of(ur))
                            mol_feat.append(fi)
                        else:
                            c["fx_incell_no_UR"] += 1

            # ---- probe-assignment universe ----
            if r.is_secondary or r.is_supplementary:
                c["skip_sec_supp"] += 1
                continue
            seq = r.query_sequence
            if not r.has_tag("pr"):
                c["no_pr"] += 1
                continue
            if r.is_unmapped:
                c["pr_but_unmapped"] += 1
                continue
            if not seq or len(seq) != 2 * HALF:
                c["no_seq50"] += 1
                continue
            pr = str(r.get_tag("pr")).split(";")[0]
            probe = probes.get(pr)
            if probe is None:
                c["probe_missing"] += 1
                continue
            gene = pr.split("|")[0]
            o = revcomp(seq) if r.is_reverse else seq
            lhs = hamming(o[:HALF], probe[:HALF])
            rhs = hamming(o[HALF:], probe[HALF:])
            prox = (lhs == 0 and rhs == 0) or any(v in gene_kmers[gene] for v in kmer_ints(o, K))
            lhs_a.append(lhs)
            rhs_a.append(rhs)
            gene_a.append(gidx[gene])
            prox_a.append(prox)
            b1.append(f"@{u}\n{encode_ordinal(base | u, g['CB_LEN'])}{r1_tail}\n+\n{q1}\n")
            b2.append(f"@{u}\n{o}\n+\n{'F' * len(o)}\n")
            u += 1
            if len(b1) >= 100_000:
                w1.write("".join(b1).encode())
                w2.write("".join(b2).encode())
                b1.clear()
                b2.clear()
    w1.write("".join(b1).encode())
    w2.write("".join(b2).encode())
    for p in gz:
        p.stdin.close()
        if p.wait() != 0:
            raise SystemExit("pigz failed")
    if u >= 1 << local_bits:
        raise SystemExit(f"{u} universe reads do not fit {local_bits} ordinal bits; use more chunks")

    c["universe"] = u
    c["cr_molecule_reads"] = len(mol_key)
    np.savez(f"{args.out}.univ.npz", lhs_sub=np.frombuffer(lhs_a, np.uint8),
             rhs_sub=np.frombuffer(rhs_a, np.uint8), gene=np.frombuffer(gene_a, np.uint32),
             proxy_keep=np.frombuffer(prox_a, np.uint8))
    np.savez(f"{args.out}.crmol.npz", key=np.frombuffer(mol_key, np.uint64),
             feature=np.frombuffer(mol_feat, np.uint32))
    out = {"bam": args.bam, "k": K, "n_records": n, "limit": args.limit,
           "chunk": args.chunk, "nchunks": args.nchunks, "local_bits": local_bits,
           "elapsed_s": round(time.time() - t0, 1), "n_called_cells": len(cells),
           "genes": genes, "features": features, "counters": dict(c)}
    with open(f"{args.out}.json", "w") as f:
        json.dump(out, f)
    print(f"[{args.bam}] {n:,} records, universe {u:,}, CR molecule reads {len(mol_key):,} "
          f"in {out['elapsed_s']} s", flush=True)


if __name__ == "__main__":
    main()
