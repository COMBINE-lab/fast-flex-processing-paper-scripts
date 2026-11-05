#!/usr/bin/env python3
"""Report how much of each file is resident in the page cache (rootless fincore via mincore(2)).

Usage: page_resident.py FILE [FILE ...]
Prints one line per file plus a TOTAL line: resident_bytes / size_bytes (pct).
Used to PROVE that lib/evict_cache.py actually drops the timed inputs before every rep.
"""
import ctypes, ctypes.util, mmap, os, sys
import numpy as np

libc = ctypes.CDLL(ctypes.util.find_library("c"), use_errno=True)
libc.mmap.restype = ctypes.c_void_p
libc.mmap.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_long]
libc.munmap.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
libc.mincore.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_ubyte)]
PAGE = os.sysconf("SC_PAGE_SIZE")
CHUNK = 1 << 30  # map 1 GiB at a time so huge FASTQs do not need one giant mapping


def resident(path):
    size = os.stat(path).st_size
    fd = os.open(path, os.O_RDONLY)
    res = 0
    try:
        off = 0
        while off < size:
            ln = min(CHUNK, size - off)
            addr = libc.mmap(None, ln, mmap.PROT_READ, mmap.MAP_SHARED, fd, off)
            if addr in (None, ctypes.c_void_p(-1).value):
                raise OSError(ctypes.get_errno(), f"mmap failed on {path}")
            npages = (ln + PAGE - 1) // PAGE
            vec = (ctypes.c_ubyte * npages)()
            try:
                if libc.mincore(addr, ln, vec) != 0:
                    raise OSError(ctypes.get_errno(), f"mincore failed on {path}")
                res += int((np.frombuffer(vec, dtype=np.uint8) & 1).sum()) * PAGE
            finally:
                libc.munmap(addr, ln)
            off += ln
    finally:
        os.close(fd)
    return min(res, size), size


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: page_resident.py FILE [FILE ...]")
    tr = ts = 0
    for f in sys.argv[1:]:
        r, s = resident(f)
        tr += r; ts += s
        print(f"{r/2**30:9.3f} / {s/2**30:9.3f} GiB  ({100*r/max(s,1):6.2f}%)  {f}")
    print(f"TOTAL {tr/2**30:.3f} / {ts/2**30:.3f} GiB resident ({100*tr/max(ts,1):.2f}%)")
