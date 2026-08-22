#!/usr/bin/env python3
"""Which libraries chunk, and which output streams their chunks are merged from.

test_orblib_chunk_merge.py checks that merging pm_hist chunks reproduces the
single-process file. It says nothing about whether a proper-motion library is
ever offered to the merger in the first place, or whether the merger is asked
for the right set of streams - which is the rest of the change that enabled
this, and is pure decision logic in orblib.py.

The two must agree. can_chunk_orbits deciding yes while chunk_kinds names a
stream the merger has no header count for would fail at merge time, after the
integration has already been spent; naming one the Fortran never wrote would
fail as a missing chunk.

Usage:
    python test_orblib_chunk_kinds.py
"""

import os
import sys

# the repo, not whatever is installed in site-packages
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dynamite import orblib_chunks as oc                   # noqa: E402
from dynamite.orblib import LegacyOrbitLibrary as OrbLib    # noqa: E402


class Log:
    """can_chunk_orbits logs its reasons; collect them instead of printing."""

    def __init__(self):
        self.messages = []

    def info(self, msg):
        self.messages.append(msg)

    debug = warning = error = info


class Stub:
    """Just enough of an orbit library for the two decision methods."""

    def __init__(self, n_hist1d=1, n_hist2d=0, n_pops=0, pops_share_aper=True,
                 legacy=False, starting_orbit=1, number_orbits=-1):
        self.settings = {'nE': 10, 'nI2': 6, 'nI3': 5, 'dithering': 2,
                         'starting_orbit': starting_orbit,
                         'number_orbits': number_orbits}
        self.n_hist1d = n_hist1d
        self.n_hist2d = n_hist2d
        self.LegacyWeightSolver = legacy
        self.mod_dir = 'stub/'
        self.logger = Log()

        # population_data entries with kin_aper None are the ones that make the
        # Fortran write a pops stream of its own, which is what still blocks
        # chunking
        class Pop:
            def __init__(self, kin_aper):
                self.kin_aper = kin_aper

        self.stars = type('Stars', (), {})()
        self.stars.population_data = [
            Pop(0 if pops_share_aper else None) for _ in range(n_pops)]

    n_orbit_starting_points = OrbLib.n_orbit_starting_points
    chunk_kinds = OrbLib.chunk_kinds
    can_chunk_orbits = OrbLib.can_chunk_orbits


def check_proper_motions_may_chunk():
    """The point of the change: PM libraries were refused outright."""
    # LOS + proper motions, the omega Cen setup
    lib = Stub(n_hist1d=3, n_hist2d=2)
    assert lib.can_chunk_orbits(), \
        'a proper-motion library is still being refused: ' \
        f'{lib.logger.messages}'
    assert lib.chunk_kinds() == ('qgrid', 'losvd_hist', 'pm_hist'), \
        lib.chunk_kinds()


def check_streams_match_what_the_fortran_writes():
    """chunk_kinds must name exactly the files write_orblib_dot_in asks for.

    write_orblib_dot_in emits the losvd file name only when n_hist1d > 0 and
    the pm file name only when n_hist2d > 0, and the Fortran creates only the
    files it was given a name for. A stream named here that was never written
    fails the merge as a missing chunk.
    """
    cases = {
        (1, 0): ('qgrid', 'losvd_hist'),               # LOS only
        (0, 1): ('qgrid', 'pm_hist'),                  # proper motions only
        (3, 2): ('qgrid', 'losvd_hist', 'pm_hist'),    # both
    }
    for (n1d, n2d), want in cases.items():
        got = Stub(n_hist1d=n1d, n_hist2d=n2d).chunk_kinds()
        assert got == want, f'n_hist1d={n1d} n_hist2d={n2d}: {got} != {want}'


def check_every_named_stream_can_be_merged():
    """Ties chunk_kinds to the merger's header table.

    These are edited in different files; a stream added to one and not the
    other is only discovered after an integration has been spent.
    """
    for n1d in (0, 1, 3):
        for n2d in (0, 1, 2):
            if n1d == 0 and n2d == 0:
                continue                    # no kinematics at all
            for kind in Stub(n_hist1d=n1d, n_hist2d=n2d).chunk_kinds():
                assert kind in oc.HEADER_RECORDS, \
                    f'{kind} has no HEADER_RECORDS entry, so it cannot merge'


def check_the_remaining_refusals_still_refuse():
    """PM was the only gate removed; the other three are real constraints."""
    # populations writing their own 0d stream: not merged
    lib = Stub(n_pops=2, pops_share_aper=False)
    assert not lib.can_chunk_orbits(), 'a populations library must not chunk'

    # populations that reuse a kinematic aperture write no pops stream
    lib = Stub(n_pops=2, pops_share_aper=True)
    assert lib.can_chunk_orbits(), \
        'populations sharing a kinematic aperture need not block chunking'

    # triaxmass reads the merged qgrid file, which does not exist yet
    lib = Stub(legacy=True)
    assert not lib.can_chunk_orbits(), 'LegacyWeightSolver must not chunk'
    assert any('triaxmass' in m for m in lib.logger.messages), \
        lib.logger.messages

    # a deliberately partial library would silently come back complete
    lib = Stub(starting_orbit=5)
    assert not lib.can_chunk_orbits(), 'a partial library must not chunk'
    lib = Stub(number_orbits=10)
    assert not lib.can_chunk_orbits(), 'a truncated library must not chunk'

    # number_orbits given explicitly but covering everything is still full
    lib = Stub(number_orbits=10 * 6 * 5)
    assert lib.can_chunk_orbits(), \
        f'an explicitly-full library must chunk: {lib.logger.messages}'


def main():
    check_proper_motions_may_chunk()
    check_streams_match_what_the_fortran_writes()
    check_every_named_stream_can_be_merged()
    check_the_remaining_refusals_still_refuse()
    print('OK')
    return 0


if __name__ == '__main__':
    sys.exit(main())
