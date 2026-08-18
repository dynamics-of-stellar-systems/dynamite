"""Check construct_nnls_matrix_and_rhs leaves the orbit library unmutated.

It zeroes the first and last velocity bin of each 1D histogram in place, and
now restores them before returning. That restore is what lets chi2_kinmap
reuse the library instead of paying for a second full read. If the restore is
exact, chi2_kinmap sees byte-identical input either way.

Also confirms the restore does not disturb the NNLS matrix itself, by
comparing against a fingerprint written by _real_orblib_check.py.

Run from ~/research/omegacen/dynamite_dataprep with the `main` conda python:
    python .../_real_chi2_kinmap_check.py [reference.npz]
"""
import sys
import time

import numpy as np

import dynamite as dyn

CONFIG = 'NGC5139_config_veldist_combined_bigger.yaml'


def main(reference=None):
    c = dyn.config_reader.Configuration(CONFIG, reset_logging=True)
    model = c.all_models.get_model_from_row(0)

    orblib = model.get_orblib()
    t0 = time.perf_counter()
    orblib.read_vel_histograms()
    t_read = time.perf_counter() - t0
    print(f'orblib read {t_read:.1f}s', flush=True)

    solver = dyn.weight_solvers.NNLS(config=c, model=model)

    before = [h.y.copy() for h in orblib.vel_histograms]
    A, b = solver.construct_nnls_matrix_and_rhs(orblib)

    ok = True
    for i, (was, hist) in enumerate(zip(before, orblib.vel_histograms)):
        same = np.array_equal(was, hist.y)
        print(f'  hist {i} {was.shape}: restored exactly = {same}')
        ok &= same
    del before

    if reference is not None:
        ref = np.load(reference)
        print(f'A unchanged vs {reference}: '
              f'{np.array_equal(ref["A"], A)} '
              f'(b: {np.array_equal(ref["b"], b)})')

    print(f'\nthe re-read chi2_kinmap would have paid for: {t_read:.1f}s')
    if not ok:
        print('FAILED: orblib was left mutated')
        sys.exit(1)
    print('OK: orblib returned unmutated, chi2_kinmap can reuse it')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
