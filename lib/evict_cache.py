#!/usr/bin/env python3
"""Rootless page-cache eviction for specific files via posix_fadvise(POSIX_FADV_DONTNEED).

Unlike `echo 3 > /proc/sys/vm/drop_caches` (global, needs root), this drops only the named
files' clean pages from the page cache, using an advisory call available to any user who can
open the file read-only. Used by scaling_study.sh to force cold-cache reads before each timed
run so the network-vs-local read comparison is honest.

Usage: evict_cache.py FILE [FILE ...]   (missing files are skipped with a warning)

Note: POSIX_FADV_DONTNEED is advisory. It reliably evicts clean, unmapped pages of a file on
local and most network filesystems; if a filesystem ignores it, the reads won't be cold (the
study prints what it evicted so you can sanity-check timing variance across reps).
"""
import os, sys

def evict(path):
    try:
        fd = os.open(path, os.O_RDONLY)
    except OSError as e:
        print(f"[evict] skip {path}: {e}", file=sys.stderr); return 0
    try:
        os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)  # offset0,len0 => whole file
        return os.fstat(fd).st_size
    finally:
        os.close(fd)

if __name__ == "__main__":
    files = sys.argv[1:]
    if not files:
        sys.exit("usage: evict_cache.py FILE [FILE ...]")
    os.sync()  # flush dirty pages first so DONTNEED can drop them
    total = sum(evict(f) for f in files)
    print(f"[evict] requested DONTNEED on {len(files)} files (~{total/1e9:.1f} GB)")
