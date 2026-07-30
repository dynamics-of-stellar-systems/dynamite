#!/usr/bin/env python3
"""A failed orbit library read must not poison the next model.

Reads are done inside a reused pathos pool worker, after ``os.chdir`` into the
model directory and after decompressing the library to a temporary file. If a
read raises, the worker used to be left in the wrong directory (so the next
model resolved ``datfil/...`` against the previous model's path) and the
decompressed library - gigabytes for a large one - stayed on disk.

Also checks that chunking is refused for the configurations it cannot produce a
correct library for.

Needs the NGC6278_out_d2 tree that test_orblib_chunking.py uses.

Usage:
    python test_orblib_failure_hygiene.py
"""

import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import dynamite as dyn                                        # noqa: E402


def make_orblib(cfg_name):
    """A LegacyOrbitLibrary for the first model of ``cfg_name``."""
    c = dyn.config_reader.Configuration(cfg_name, reset_logging=True,
                                        reset_existing_output=False)
    mod = c.all_models.get_model_from_row(0)
    return dyn.orblib.LegacyOrbitLibrary(config=c, mod_dir=mod.directory_noml,
                                         parset=mod.parset)


def check_failed_read_is_clean(lib):
    """A read that cannot decompress must restore cwd and leave no temp file."""
    before = os.getcwd()
    datfil = os.path.join(lib.mod_dir, 'datfil')
    existing = set(glob.glob(os.path.join(datfil, '*.dat')))

    try:
        # no such family, so bunzip2 has nothing to read and _decompress raises
        lib.read_orbit_base('nosuchorblib')
    except RuntimeError as e:
        assert 'cannot decompress' in str(e), str(e)
    else:
        raise AssertionError('a missing orbit library was read without error')

    assert os.getcwd() == before, \
        f'working directory left at {os.getcwd()}, expected {before}'
    leaked = set(glob.glob(os.path.join(datfil, '*.dat'))) - existing
    assert not leaked, f'decompressed files left behind: {sorted(leaked)}'


def check_real_read_still_works(lib):
    """The decorators must not disturb a successful read."""
    before = os.getcwd()
    datfil = os.path.join(lib.mod_dir, 'datfil')
    existing = set(glob.glob(os.path.join(datfil, '*.dat')))
    _, density_3d = lib.read_orbit_base('orblib')
    n_orb = (lib.settings['nE'] * lib.settings['nI2'] * lib.settings['nI3'])
    assert density_3d.shape[0] == n_orb, \
        f'read {density_3d.shape[0]} orbits, expected {n_orb}'
    assert os.getcwd() == before, 'working directory not restored'
    leaked = set(glob.glob(os.path.join(datfil, '*.dat'))) - existing
    assert not leaked, f'decompressed files left behind: {sorted(leaked)}'


def check_chunking_refused_when_unsafe(lib):
    """can_chunk_orbits must refuse what the chunked path cannot do."""
    assert lib.can_chunk_orbits(), \
        'this configuration should be chunkable; the rest of the check is moot'

    # triaxmass reads the merged datfil/orblib_qgrid.dat, which does not exist
    # until after the chunked script has already finished
    lib.LegacyWeightSolver = True
    assert not lib.can_chunk_orbits(), 'chunked LegacyWeightSolver not refused'
    lib.LegacyWeightSolver = False

    # the chunks' ranges replace these outright, so a deliberately partial
    # library would silently come back complete
    for key, bad in (('starting_orbit', 5), ('number_orbits', 10)):
        keep = lib.settings[key]
        lib.settings[key] = bad
        assert not lib.can_chunk_orbits(), f'chunked partial library ({key}={bad}) not refused'
        lib.settings[key] = keep

    assert lib.can_chunk_orbits(), 'checks did not restore the settings'


def check_stale_ics_are_detected(lib):
    """begin.dat for a different orbit grid must not be reused.

    Model directories are orblib_<iteration>_<row>, so changing nE/nI2/nI3/
    dithering and re-running into the same output directory finds initial
    conditions for the old grid. The Fortran then either refuses them with
    "Not so many orbits" or integrates the wrong starting points.
    """
    assert lib.ics_match_settings(), \
        'the seeded tree should match its own config; the rest is moot'
    for key in ('nE', 'nI2', 'nI3', 'dithering'):
        keep = lib.settings[key]
        lib.settings[key] = keep + 1
        assert not lib.ics_match_settings(), f'stale ICS not detected for {key}'
        lib.settings[key] = keep
    assert lib.ics_match_settings(), 'checks did not restore the settings'


def check_truncated_losvd_is_named(lib):
    """A short losvd file must fail with a message naming it, not IndexError."""
    import numpy as np
    datfil = os.path.join(lib.mod_dir, 'datfil')
    src = os.path.join(datfil, 'orblib_losvd_hist.dat.bz2')
    trunc = os.path.join(datfil, 'truncated_losvd.dat')
    # decompress, then keep only the header and a couple of pairs
    lib._decompressed = []
    lib._decompress(src, trunc)
    buf = np.fromfile(trunc, dtype=np.uint8)
    open(trunc, 'wb').write(buf[:200].tobytes())

    n_ap = 3
    try:
        lib._read_losvd_hist_vectorised(
            trunc, 10 ** 6, np.zeros(n_ap, dtype=int), np.zeros(1, dtype=int),
            [21], [np.zeros((10 ** 6, 21, n_ap))])
    except ValueError as e:
        assert 'truncated' in str(e) and 'truncated_losvd.dat' in str(e), str(e)
    except IndexError as e:
        raise AssertionError(f'bare IndexError instead of a named error: {e}')
    else:
        raise AssertionError('a truncated losvd file was read without error')
    finally:
        os.remove(trunc)


def main():
    os.chdir(HERE)
    lib = make_orblib('user_test_config_ml_dither2.yaml')
    check_chunking_refused_when_unsafe(lib)
    check_stale_ics_are_detected(lib)
    check_failed_read_is_clean(lib)
    check_truncated_losvd_is_named(lib)
    check_real_read_still_works(lib)
    print('OK')
    return 0


if __name__ == '__main__':
    sys.exit(main())
