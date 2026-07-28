#!/usr/bin/env python3
"""Merge per-chunk orbit library files into the single files dynamite expects.

Both output streams are sequences of Fortran unformatted records: a fixed
header, then a per-orbit (or per-orbit-per-aperture) body. Chunks cover
disjoint, consecutive orbit ranges, so merging is: keep the first chunk's
header, then concatenate every chunk's body in order.

    _qgrid.dat       header = 5 records, and record 0 carries the orbit count
    _losvd_hist.dat  header = 1 record, no orbit count to patch
"""
import os
import struct
import sys

import numpy as np

HEADER_RECORDS = {"qgrid": 5, "losvd_hist": 1}


def record_offsets(buf, n_records):
    """Byte offset just past the first n_records Fortran records."""
    p = 0
    for _ in range(n_records):
        if p + 4 > len(buf):
            raise ValueError("truncated file while reading header")
        n = int(np.frombuffer(buf[p:p + 4], dtype=np.int32)[0])
        trail = int(np.frombuffer(buf[p + 4 + n:p + 8 + n], dtype=np.int32)[0])
        if trail != n:
            raise ValueError(f"record markers disagree: {n} vs {trail}")
        p += 8 + n
    return p


def footer_start(buf):
    """Offset of the final record, which every file ends with.

    The integrator closes each file with a 1-byte record, so a merged file
    must carry exactly one of them rather than one per chunk. Located from the
    trailing length marker, so this costs nothing on a large file.
    """
    n = int(np.frombuffer(buf[-4:], dtype=np.int32)[0])
    start = len(buf) - 8 - n
    lead = int(np.frombuffer(buf[start:start + 4], dtype=np.int32)[0])
    if lead != n:
        raise ValueError(f"trailing record markers disagree: {lead} vs {n}")
    return start


def merge(chunk_paths, out_path, kind, total_orbits=None):
    n_hdr = HEADER_RECORDS[kind]
    bufs = [np.fromfile(p, dtype=np.uint8) for p in chunk_paths]
    head_end = record_offsets(bufs[0], n_hdr)
    header = bytearray(bufs[0][:head_end].tobytes())

    if kind == "qgrid" and total_orbits is not None:
        # first record's payload starts at byte 4; its first int32 is norb
        struct.pack_into("<i", header, 4, total_orbits)

    with open(out_path, "wb") as f:
        f.write(header)
        for buf in bufs:
            f.write(buf[record_offsets(buf, n_hdr):footer_start(buf)].tobytes())
        f.write(bufs[-1][footer_start(bufs[-1]):].tobytes())
    return out_path


def merge_model(datfil, base, n_chunks, total_orbits):
    """Merge {base}_c{k}_* chunk files into {base}_*."""
    for kind, suffix in (("qgrid", "_qgrid.dat"),
                         ("losvd_hist", "_losvd_hist.dat")):
        parts = [os.path.join(datfil, f"{base}_c{k}{suffix}")
                 for k in range(n_chunks)]
        missing = [p for p in parts if not os.path.exists(p)]
        if missing:
            raise FileNotFoundError(missing[0])
        merge(parts, os.path.join(datfil, base + suffix), kind, total_orbits)


if __name__ == "__main__":
    datfil, base, n, norb = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    merge_model(datfil, base, n, norb)
    print(f"merged {n} chunks -> {base}_qgrid.dat, {base}_losvd_hist.dat")
