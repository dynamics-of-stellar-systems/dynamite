#!/usr/bin/env python3
"""Chunk ranges must cover the whole library, at any dithering.

The end-to-end chunking test (test_orblib_chunking.py) runs with dithering: 1,
where an off-by-dithering^3 error in the starting-point count is invisible.
This check is the cheap one that isn't.

Usage:
    python test_orbit_chunk_bounds.py
"""

import os
import sys

# the repo, not whatever is installed in site-packages
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dynamite.orblib import LegacyOrbitLibrary as OrbLib     # noqa: E402


class Stub:
    """Just enough of an orbit library for the two arithmetic methods."""

    def __init__(self, nE, nI2, nI3, dithering):
        self.settings = {'nE': nE, 'nI2': nI2, 'nI3': nI3,
                         'dithering': dithering}

    n_orbit_starting_points = OrbLib.n_orbit_starting_points
    orbit_chunk_bounds = OrbLib.orbit_chunk_bounds


def main():
    # the omega Cen configuration that exposed this, plus the dithering: 1 and
    # does-not-divide-evenly cases
    for nE, nI2, nI3, dith in ((40, 20, 15, 3), (10, 6, 5, 1), (7, 5, 3, 2)):
        lib = Stub(nE, nI2, nI3, dith)
        # read_orbit_base expects exactly this many orbits in the merged file
        expected = nE * nI2 * nI3
        assert lib.n_orbit_starting_points() == expected, \
            f'{nE}x{nI2}x{nI3} d{dith}: {lib.n_orbit_starting_points()}'
        for n_chunks in (1, 2, 5, 7, expected, expected + 3):
            bounds = lib.orbit_chunk_bounds(n_chunks)
            assert sum(n for _, n, _ in bounds) == expected, \
                f'chunks do not cover the library: {n_chunks} chunks'
            # consecutive and 1-based, as the Fortran input file wants
            start = 1
            for s, n, _ in bounds:
                assert s == start, f'gap or overlap at chunk starting {s}'
                start += n
    print('OK')
    return 0


if __name__ == '__main__':
    sys.exit(main())
