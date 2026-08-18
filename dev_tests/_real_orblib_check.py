"""Build the NNLS matrix from a REAL orbit library, profile it, and dump a
fingerprint so two git revisions can be compared.

Run from ~/research/omegacen/dynamite_dataprep with the `main` conda python:
    python .../_real_orblib_check.py <out.npz>
"""
import cProfile
import pstats
import sys
import time

import numpy as np

import dynamite as dyn

CONFIG = 'NGC5139_config_veldist_combined_bigger.yaml'


def main(outfile):
    t0 = time.perf_counter()
    c = dyn.config_reader.Configuration(CONFIG, reset_logging=True)
    mods = c.all_models
    model = mods.get_model_from_row(0)
    print(f'config+model loaded in {time.perf_counter()-t0:.1f}s', flush=True)

    t0 = time.perf_counter()
    orblib = model.get_orblib()
    orblib.read_vel_histograms()
    print(f'orblib read in {time.perf_counter()-t0:.1f}s', flush=True)
    print('n_orbs', orblib.n_orbs,
          'vel_hist shapes', [h.y.shape for h in orblib.vel_histograms],
          flush=True)

    solver = dyn.weight_solvers.NNLS(config=c, model=model)

    pr = cProfile.Profile()
    t0 = time.perf_counter()
    pr.enable()
    A, b = solver.construct_nnls_matrix_and_rhs(orblib)
    pr.disable()
    wall = time.perf_counter() - t0
    print(f'construct_nnls_matrix_and_rhs: {wall:.2f}s  A={A.shape} {A.dtype}',
          flush=True)
    pstats.Stats(pr).sort_stats('cumulative').print_stats(18)

    np.savez(outfile, A=A, b=b, wall=wall)
    print('A checksum', float(np.abs(A).sum()), 'b checksum', float(np.abs(b).sum()))


if __name__ == '__main__':
    main(sys.argv[1])
