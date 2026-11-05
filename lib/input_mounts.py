#!/usr/bin/env python3
"""Fail-loud check that every input of a timed run is on the local NVMe scratch filesystem.

Every PATH is resolved with realpath (a directory: every file under it, symlinks followed) and its
filesystem (st_dev) must equal the one of --scr. The table path -> realpath -> device -> mount source
(findmnt -T) is written to --out (TSV); one summary line goes to stdout:
    "<n> input files on <source> (<target>)"
Exit 1 lists every input that is missing or on another filesystem (for example network storage or NFS).

Usage: input_mounts.py --scr DIR --out TSV PATH...
"""
import argparse
import os
import subprocess
import sys


def files_of(p):
    if os.path.isdir(p):
        for root, _, names in os.walk(p, followlinks=True):
            for n in names:
                yield os.path.join(root, n)
    else:
        yield p


def mount_of(p):
    o = subprocess.check_output(["findmnt", "-n", "-o", "SOURCE,TARGET", "-T", p], text=True).split()
    return o[0], o[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scr", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("paths", nargs="+")
    a = ap.parse_args()
    want = os.stat(a.scr).st_dev
    scr_src, scr_tgt = mount_of(a.scr)
    rows, bad, srcs = [], [], {}
    for p in a.paths:
        if not os.path.exists(p):
            bad.append(f"{p}: missing")
            continue
        for f in files_of(p):
            r = os.path.realpath(f)
            if not os.path.exists(r):
                bad.append(f"{f} -> {r}: missing")
                continue
            dev = os.stat(r).st_dev
            if dev not in srcs:
                srcs[dev] = mount_of(r)
            src, tgt = srcs[dev]
            rows.append((f, r, f"{os.major(dev)}:{os.minor(dev)}", src, tgt))
            if dev != want:
                bad.append(f"{f} -> {r}: on {src} ({tgt}), not on {scr_src} ({scr_tgt}) like {a.scr}")
    with open(a.out, "w") as fh:
        fh.write("path\trealpath\tdevice\tsource\ttarget\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
    if bad:
        for b in bad:
            print(f"input_mounts !! {b}", file=sys.stderr)
        sys.exit(1)
    print(f"{len(rows)} input files on {scr_src} ({scr_tgt})")


if __name__ == "__main__":
    main()
