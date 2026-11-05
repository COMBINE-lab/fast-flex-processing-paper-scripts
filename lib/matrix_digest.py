#!/usr/bin/env python3
"""Order-independent content digest of an alevin-fry count matrix (quants_mat.mtx[.gz] + rows/cols).

The digest is a sum (mod 2^64) of a 64-bit mix of (row-name hash, col-name hash, value bits)
over every stored entry, so it is invariant to row/column ORDER and to entry order: two
matrices that hold the same (barcode, gene, count) triples get the same digest, whatever
thread count or decoder produced them. Explicit zeros are skipped.

Usage: matrix_digest.py ALEVIN_DIR   (dir holding quants_mat.mtx[.gz], quants_mat_rows.txt, quants_mat_cols.txt)
Prints one tab-separated line: n_rows n_cols nnz total digest
"""
import gzip, hashlib, os, sys
import numpy as np
import pandas as pd

M64 = np.uint64(0xFFFFFFFFFFFFFFFF)


def name_hashes(path):
    with open(path) as fh:
        names = [l.rstrip("\n") for l in fh]
    return np.array([int.from_bytes(hashlib.blake2b(n.encode(), digest_size=8).digest(), "little")
                     for n in names], dtype=np.uint64)


def mix(x):  # splitmix64 finalizer; uint64 arithmetic wraps
    x = x ^ (x >> np.uint64(30)); x = x * np.uint64(0xBF58476D1CE4E5B9)
    x = x ^ (x >> np.uint64(27)); x = x * np.uint64(0x94D049BB133111EB)
    return x ^ (x >> np.uint64(31))


def main(d):
    mtx = next(os.path.join(d, n) for n in ("quants_mat.mtx", "quants_mat.mtx.gz") if os.path.exists(os.path.join(d, n)))
    rh = name_hashes(os.path.join(d, "quants_mat_rows.txt"))
    ch = name_hashes(os.path.join(d, "quants_mat_cols.txt"))
    op = gzip.open if mtx.endswith(".gz") else open
    with op(mtx, "rt") as fh:
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        nr, nc, nnz_hdr = (int(x) for x in line.split())
        if (nr, nc) != (len(rh), len(ch)):
            sys.exit(f"shape {nr}x{nc} != rows/cols files {len(rh)}x{len(ch)}")
        acc = np.uint64(0); nnz = 0; total = 0.0
        with np.errstate(over="ignore"):
            for chunk in pd.read_csv(fh, sep=r"\s+", header=None, names=["i", "j", "v"],
                                     dtype={"i": np.int64, "j": np.int64, "v": np.float64},
                                     chunksize=20_000_000, engine="c"):
                v = chunk["v"].to_numpy()
                keep = v != 0
                i = chunk["i"].to_numpy()[keep] - 1; j = chunk["j"].to_numpy()[keep] - 1; v = v[keep]
                h = mix(rh[i] ^ mix(ch[j] + np.uint64(0x9E3779B97F4A7C15)) ^ mix(v.view(np.uint64)))
                acc = np.uint64((int(acc) + int(h.sum(dtype=np.uint64))) & 0xFFFFFFFFFFFFFFFF)
                nnz += int(keep.sum()); total += float(v.sum())
    tot = int(total) if total == int(total) else total
    print(f"{nr}\t{nc}\t{nnz}\t{tot}\t{int(acc):016x}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
