#!/usr/bin/env python3
"""Run an orbit library monolithically and in N chunks, and compare.

Drives legacy_fortran/orblib_new_mirror directly by rewriting the orbit-range
and output-filename lines of an existing orblib.in.

    python chunk_test.py <model_dir> --chunks 4 [--tag T] [--norb N] [--parallel]
"""
import argparse
import os
import re
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, "/Users/pesmith/research/dynamite")
from dynamite.orblib import LegacyOrbitLibrary as L        # noqa: E402

EXE = os.path.expanduser("~/research/dynamite/legacy_fortran/orblib_new_mirror")


class _Parser:
    """Borrow the vectorised record parser without constructing an orblib."""
    _walk_losvd_records = L._walk_losvd_records
    _gather_float64 = staticmethod(L._gather_float64)
    _read_losvd_hist_vectorised = L._read_losvd_hist_vectorised


def write_in(src, dst, start, number, tag):
    """Copy an orblib.in, setting the orbit range and output filenames."""
    lines = open(src).read().splitlines()
    out = []
    for line in lines:
        if re.search(r"\[starting orbit\]", line):
            line = re.sub(r"^\s*\d+", str(start), line, count=1)
        elif re.search(r"orbits\s+to\s+int", line):
            line = re.sub(r"^\s*-?\d+", str(number), line, count=1)
        elif "_qgrid.dat" in line and line.strip().startswith('"'):
            line = re.sub(r'"[^"]*_qgrid\.dat"', f'"datfil/{tag}_qgrid.dat"', line)
        elif "_losvd_hist.dat" in line and line.strip().startswith('"'):
            line = re.sub(r'"[^"]*_losvd_hist\.dat"',
                          f'"datfil/{tag}_losvd.dat"', line)
        elif "orbclass.out" in line and line.strip().startswith('"'):
            line = re.sub(r'"[^"]*orbclass\.out"', f'"datfil/{tag}_class.out"',
                          line)
        out.append(line)
    open(dst, "w").write("\n".join(out) + "\n")


def run(model_dir, tag, start, number, src_in, background=False):
    for suffix in ("_qgrid.dat", "_qgrid.dat.tmp", "_losvd.dat", "_class.out"):
        p = os.path.join(model_dir, "datfil", tag + suffix)
        if os.path.exists(p):
            os.remove(p)
    in_path = os.path.join(model_dir, "infil", tag + ".in")
    write_in(src_in, in_path, start, number, tag)
    log = open(os.path.join(model_dir, "datfil", tag + ".log"), "w")
    proc = subprocess.Popen([EXE], stdin=open(in_path), stdout=log,
                            stderr=subprocess.STDOUT, cwd=model_dir)
    if background:
        return proc
    proc.wait()
    return proc


def hdr(path):
    b = np.memmap(path, dtype=np.uint8, mode="r")
    i = b[:16].view(np.int32)
    return int(i[1]), int(i[2])          # n_apertures, nvhist


def read_losvd(path, norb):
    n_ap, nvh = hdr(path)
    nv = 2 * nvh + 1
    out = [np.zeros((norb, nv, n_ap))]
    _Parser()._read_losvd_hist_vectorised(
        path, norb, np.zeros(n_ap, int), np.array([0]), [nv], out)
    return out[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("model_dir")
    ap.add_argument("--chunks", type=int, default=4)
    ap.add_argument("--norb", type=int, required=True)
    ap.add_argument("--tag", default="CT")
    ap.add_argument("--parallel", action="store_true")
    ap.add_argument("--src", default="infil/orblib.in")
    args = ap.parse_args()

    md = os.path.abspath(os.path.expanduser(args.model_dir))
    src = os.path.join(md, args.src)
    N, norb = args.chunks, args.norb

    # monolithic
    t0 = time.perf_counter()
    run(md, f"{args.tag}mono", 1, -1, src)
    t_mono = time.perf_counter() - t0
    print(f"monolithic: {t_mono:.1f}s", flush=True)

    # chunks
    base, rem = divmod(norb, N)
    bounds, s = [], 1
    for k in range(N):
        n = base + (1 if k < rem else 0)
        bounds.append((s, n))
        s += n
    t0 = time.perf_counter()
    if args.parallel:
        procs = [run(md, f"{args.tag}c{k}", st, n, src, background=True)
                 for k, (st, n) in enumerate(bounds)]
        for p in procs:
            p.wait()
    else:
        for k, (st, n) in enumerate(bounds):
            run(md, f"{args.tag}c{k}", st, n, src)
    t_chunk = time.perf_counter() - t0
    print(f"{N} chunks{' (parallel)' if args.parallel else ''}: "
          f"{t_chunk:.1f}s  -> {t_mono / max(t_chunk, 1e-9):.2f}x", flush=True)

    # compare
    mono = read_losvd(os.path.join(md, "datfil", f"{args.tag}mono_losvd.dat"),
                      norb)
    ok = True
    off = 0
    for k, (st, n) in enumerate(bounds):
        f = os.path.join(md, "datfil", f"{args.tag}c{k}_losvd.dat")
        got = read_losvd(f, n)
        ref = mono[off:off + n]
        same = np.array_equal(ref, got)
        ok &= same
        if same:
            print(f"  chunk {k} (orbits {st}-{st + n - 1}): identical")
        else:
            d = np.abs(ref - got).sum() / max(np.abs(ref).sum(), 1e-300)
            print(f"  chunk {k} (orbits {st}-{st + n - 1}): DIFFER relL1={d:.3e}")
        off += n
    print(f"\nALL CHUNKS IDENTICAL: {ok}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
