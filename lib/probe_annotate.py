#!/usr/bin/env python3
"""
Build the probe annotation table that every probe-level analysis depends on.

For each probe in a 10x Flex probe-set CSV, locate its 50-mer against the matching
genome + GTF and record what the probe actually measures:

  genomic_class   within_exon | junction | intronic | not_found
                  Re-derived, not taken from the CSV's `region` column: 11 % of
                  `unspliced`-labelled probes span a junction under GENCODE 2024-A
                 .
  n_tx_hit        annotated isoforms of the gene whose spliced sequence contains the probe
  n_tx_gene       annotated isoforms of the gene
  constitutive    n_tx_hit == n_tx_gene, i.e. the probe cannot distinguish isoforms
  tx_set          the isoforms hit, sorted and comma-joined (the isoform-set key)

Chromosomes are read one at a time and held in memory so the genome FASTA is streamed
sequentially — random seeks over a network (NFS) mount are ~100x slower.

Usage:
  probe_annotate.py --probe-csv CSV --gtf genes.gtf.gz --fasta genome.fa --out TSV
"""
import argparse
import gzip
from collections import defaultdict

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def read_probes(path):
    probes = defaultdict(list)
    with open(path) as f:
        hdr = None
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split(",")
            if hdr is None:
                hdr = {c: i for i, c in enumerate(p)}
                continue
            probes[p[hdr["gene_id"]]].append((
                p[hdr["probe_id"]], p[hdr["probe_seq"]],
                p[hdr["region"]], p[hdr["included"]], p[hdr["gene_name"]]))
    return probes


def read_exons(gtf, wanted):
    """gene_id -> {transcript_id: [(start, end)]}, plus (chrom, strand) per gene."""
    exons = defaultdict(lambda: defaultdict(list))
    meta = {}
    op = gzip.open if gtf.endswith(".gz") else open
    with op(gtf, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.split("\t", 9)
            if c[2] != "exon":
                continue
            a = c[8]
            i = a.find('gene_id "')
            g = a[i + 9:a.find('"', i + 9)]
            if g not in wanted:
                continue
            j = a.find('transcript_id "')
            tx = a[j + 15:a.find('"', j + 15)]
            exons[g][tx].append((int(c[3]), int(c[4])))
            meta[g] = (c[0], c[6])
    return exons, meta


def chrom_reader(fasta):
    """Yield (name, sequence) streaming the FASTA once."""
    name, buf = None, []
    with open(fasta) as f:
        for line in f:
            if line[0] == ">":
                if name:
                    yield name, "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.rstrip())
    if name:
        yield name, "".join(buf).upper()


def classify(seq, span, span_lo, gene_exons):
    """(genomic_class, tx_set) for one probe against one gene."""
    sq, rq = seq, rc(seq)
    k = span.find(sq)
    if k < 0:
        k = span.find(rq)
    if k >= 0:
        s0, e0 = span_lo + k, span_lo + k + len(seq) - 1
        in_exon = any(s0 >= s and e0 <= e
                      for parts in gene_exons.values() for s, e in parts)
        return ("within_exon" if in_exon else "intronic"), None
    return "junction", None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe-csv", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-span", type=int, default=3_000_000,
                    help="skip genes whose annotated span exceeds this (default 3 Mb)")
    args = ap.parse_args()

    probes = read_probes(args.probe_csv)
    print(f"[probe_annotate] {len(probes)} genes, "
          f"{sum(len(v) for v in probes.values())} probes", flush=True)

    exons, meta = read_exons(args.gtf, set(probes))
    print(f"[probe_annotate] exons for {len(exons)} genes", flush=True)

    by_chrom = defaultdict(list)
    for g in exons:
        by_chrom[meta[g][0]].append(g)

    out = open(args.out, "w")
    out.write("probe_id\tgene_id\tgene_name\tregion_csv\tincluded\tgenomic_class"
              "\tn_tx_hit\tn_tx_gene\tconstitutive\ttx_set\n")
    n_done = 0
    for chrom, seq in chrom_reader(args.fasta):
        genes = by_chrom.get(chrom)
        if not genes:
            continue
        for g in genes:
            ge = exons[g]
            strand = meta[g][1]
            allex = [iv for parts in ge.values() for iv in parts]
            lo, hi = min(s for s, _ in allex), max(e for _, e in allex)
            if hi - lo > args.max_span:
                span = None
            else:
                span = seq[lo - 1:hi]
            txseq = {}
            if span is not None:
                for tx, parts in ge.items():
                    parts = sorted(parts)
                    s = "".join(span[a - lo:b - lo + 1] for a, b in parts)
                    txseq[tx] = rc(s) if strand == "-" else s
            ntx = len(ge)
            for pid, pseq, region, inc, gname in probes[g]:
                if span is None:
                    cls, hits = "not_found", []
                else:
                    cls, _ = classify(pseq, span, lo, ge)
                    rq = rc(pseq)
                    hits = sorted(tx for tx, s in txseq.items()
                                  if pseq in s or rq in s)
                    if cls == "junction" and not hits:
                        cls = "not_found"
                out.write(f"{pid}\t{g}\t{gname}\t{region}\t{inc}\t{cls}\t{len(hits)}"
                          f"\t{ntx}\t{int(len(hits) == ntx and hits != [])}"
                          f"\t{','.join(hits)}\n")
            n_done += 1
        print(f"[probe_annotate] {chrom}: {len(genes)} genes ({n_done} total)", flush=True)
    out.close()
    print(f"[probe_annotate] wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
