#!/usr/bin/env python3
"""Merging chunk files must reproduce the single-process file byte for byte.

Builds synthetic orbit library files with the record layout the Fortran writes
(orblib_f_new_mirror.f90: integrator_setup_write + qgrid_setup_write, then
integrator_write/qgrid_write per orbit, then the one-byte record from
output_close), splits them into chunks, merges, and compares.

Fast and self-contained - no Fortran, no model tree. Complements
test_orblib_chunking.py, which does the same end to end but takes hours.

Usage:
    python test_orblib_chunk_merge.py
"""

import os
import struct
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dynamite import orblib_chunks as oc                      # noqa: E402


def rec(payload):
    """One Fortran unformatted record: length marker, payload, length marker."""
    n = len(payload)
    return struct.pack('<i', n) + payload + struct.pack('<i', n)


def qgrid_body(orbit, nq=4):
    """The three records integrator_write/qgrid_write emit per orbit."""
    out = rec(struct.pack('<5i', orbit, 1, 2, 3, 0))          # orbit, E1, I2, I3
    out += rec(struct.pack('<8i', *range(orbit, orbit + 8)))   # orbittypes
    out += rec(struct.pack(f'<{nq}d', *(float(orbit) + i for i in range(nq))))
    return out


def qgrid_header(n_orbits, ndith=3):
    """integrator_setup_write (1 record) + qgrid_setup_write (4 records)."""
    out = rec(struct.pack('<5i', n_orbits, 4, 4, 4, ndith))
    out += rec(struct.pack('<4i', 16, 2, 2, 2))
    for _ in range(3):
        out += rec(struct.pack('<3d', 0.0, 1.0, 2.0))
    return out


def losvd_header():
    """histogram_setup_write: one record of (nconstr, nvcube, dvcube)."""
    return rec(struct.pack('<iif', 3, 66, 1.5))


def losvd_body(orbit, n_ap=3):
    """One (ivmin, ivmax) record per aperture, plus a payload when non-empty."""
    out = b''
    for ap in range(n_ap):
        if (orbit + ap) % 3 == 0:            # empty pair: no payload record
            out += rec(struct.pack('<2i', 1, 0))
        else:
            nv = 1 + (orbit + ap) % 4
            out += rec(struct.pack('<2i', -nv // 2, -nv // 2 + nv - 1))
            out += rec(struct.pack(f'<{nv}d', *(float(orbit * 10 + i)
                                                for i in range(nv))))
    return out


CLOSING = rec(b' ')     # output_close writes a single blank


def write_family(datfil, fileroot, tag, orbits, n_orbits_header):
    """Write one chunk's (or the whole library's) three output files."""
    with open(os.path.join(datfil, f'{fileroot}{tag}_qgrid.dat'), 'wb') as f:
        f.write(qgrid_header(n_orbits_header))
        for o in orbits:
            f.write(qgrid_body(o))
        f.write(CLOSING)
    with open(os.path.join(datfil,
                           f'{fileroot}{tag}_losvd_hist.dat'), 'wb') as f:
        f.write(losvd_header())
        for o in orbits:
            f.write(losvd_body(o))
        f.write(CLOSING)
    with open(os.path.join(datfil,
                           f'{fileroot}{tag}.dat_orbclass.out'), 'w') as f:
        for o in orbits:
            f.write(f'  {o}.0  {o}.1  {o}.2\n')


def check_split(datfil, n_orbits, chunk_sizes):
    """Merge a split of ``n_orbits`` and require byte-identity with one file."""
    orbits = list(range(1, n_orbits + 1))
    write_family(datfil, 'reference', '', orbits, n_orbits)

    tags, lo = [], 0
    for k, size in enumerate(chunk_sizes):
        tag = f'_c{k}'
        tags.append(tag)
        # every chunk writes the full-library header, as the Fortran does
        write_family(datfil, 'orblib', tag, orbits[lo:lo + size], n_orbits)
        lo += size
    assert lo == n_orbits, 'test bug: chunk sizes do not sum to n_orbits'

    oc.merge_chunks(datfil, 'orblib', tags, n_orbits)

    for kind in ('qgrid', 'losvd_hist'):
        got = open(os.path.join(datfil, f'orblib_{kind}.dat'), 'rb').read()
        want = open(os.path.join(datfil, f'reference_{kind}.dat'), 'rb').read()
        assert got == want, (
            f'{kind}: merged {len(got)} bytes != reference {len(want)} bytes '
            f'for split {chunk_sizes}')
    got = open(os.path.join(datfil, 'orblib.dat_orbclass.out')).read()
    want = open(os.path.join(datfil, 'reference.dat_orbclass.out')).read()
    assert got == want, f'orbclass mismatch for split {chunk_sizes}'

    # the chunks must be gone, so a later run cannot merge a stale one
    for tag in tags:
        for kind in ('qgrid', 'losvd_hist'):
            leftover = os.path.join(datfil, f'orblib{tag}_{kind}.dat')
            assert not os.path.isfile(leftover), f'chunk left behind: {leftover}'
        leftover = os.path.join(datfil, f'orblib{tag}.dat_orbclass.out')
        assert not os.path.isfile(leftover), f'chunk left behind: {leftover}'


def check_reader_walks_merged_losvd(datfil, n_orbits, n_ap=3):
    """The vectorised reader must parse the merged losvd file.

    Ties the merge's HEADER_RECORDS['losvd_hist'] = 1 to the reader's own
    header skip: if either drifts, every payload offset is wrong. The
    investigation of the omega Cen failure suspected this header of being
    malformed because it holds (nconstr, nvcube, dvcube) rather than an orbit
    count - it does not, and this pins that down.
    """
    import numpy as np
    from dynamite.orblib import LegacyOrbitLibrary

    buf = np.fromfile(os.path.join(datfil, 'orblib_losvd_hist.dat'),
                      dtype=np.uint8)
    start, ivmin, nv = LegacyOrbitLibrary._walk_losvd_records(
        None, buf, n_orbits * n_ap)
    for orbit in range(1, n_orbits + 1):
        for ap in range(n_ap):
            k = (orbit - 1) * n_ap + ap
            if (orbit + ap) % 3 == 0:
                assert nv[k] == 0, f'orbit {orbit} ap {ap}: expected empty'
                continue
            want_nv = 1 + (orbit + ap) % 4
            assert nv[k] == want_nv, f'orbit {orbit} ap {ap}: nv {nv[k]}'
            assert ivmin[k] == -want_nv // 2, \
                f'orbit {orbit} ap {ap}: ivmin {ivmin[k]}'
            # start is an int32 index; the payload is the values we wrote
            vals = buf[start[k] * 4:start[k] * 4 + 8 * want_nv].view(np.float64)
            want = [float(orbit * 10 + i) for i in range(want_nv)]
            assert list(vals) == want, f'orbit {orbit} ap {ap}: {vals} != {want}'


def check_short_library_is_caught(datfil):
    """Chunks that do not cover the library must raise, not patch the header.

    This is the omega Cen failure: the chunk ranges covered a
    dithering^3 fraction of the library, and overwriting the header count made
    the reader run off the end of the body instead of reporting the shortfall.
    """
    # a header that is wrong outright
    write_family(datfil, 'short', '_c0', [1, 2, 3], 3)
    try:
        oc.merge_chunks(datfil, 'short', ['_c0'], 81)
    except ValueError as e:
        assert '3 orbits' in str(e) and '81' in str(e), str(e)
    else:
        raise AssertionError('a wrong header was merged without complaint')

    # and the case that actually happened: the Fortran writes the header from
    # begin.dat, so it states the full library size even when the process was
    # asked for a fraction of it. Only the body is short.
    write_family(datfil, 'shortbody', '_c0', [1, 2, 3], 81)
    try:
        oc.merge_chunks(datfil, 'shortbody', ['_c0'], 81)
    except ValueError as e:
        assert 'hold 3 orbits' in str(e) and '81' in str(e), str(e)
        return
    raise AssertionError('a short body was merged without complaint')


def main():
    with tempfile.TemporaryDirectory() as datfil:
        # evenly divisible, one chunk, and remainder spread over the chunks
        for n_orbits, sizes in ((12, [3, 3, 3, 3]),
                                (12, [12]),
                                (12, [3, 3, 3, 2, 1]),
                                (7, [2, 2, 1, 1, 1]),
                                (5, [1, 1, 1, 1, 1])):
            check_split(datfil, n_orbits, sizes)
            check_reader_walks_merged_losvd(datfil, n_orbits)
        check_short_library_is_caught(datfil)
    print('OK')
    return 0


if __name__ == '__main__':
    sys.exit(main())
