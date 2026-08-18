"""Profile a full NNLS.solve() on a real orbit library.

Run from ~/research/omegacen/dynamite_dataprep with the `main` conda python:
    python .../_real_solve_profile.py [adelie|scipy]
"""
import cProfile
import pstats
import resource
import sys
import time

import numpy as np

import dynamite as dyn

CONFIG = 'NGC5139_config_veldist_combined_bigger.yaml'


def rss_gib():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 2**30 if sys.platform == 'darwin' else r / 2**20


def main(which):
    c = dyn.config_reader.Configuration(CONFIG, reset_logging=True)
    model = c.all_models.get_model_from_row(0)

    t0 = time.perf_counter()
    orblib = model.get_orblib()
    orblib.read_vel_histograms()
    print(f'orblib read {time.perf_counter()-t0:.1f}s, peak RSS {rss_gib():.1f} GiB',
          flush=True)

    solver = dyn.weight_solvers.NNLS(config=c, model=model)
    solver.nnls_solver = which
    print('solver', solver.nnls_solver, flush=True)

    pr = cProfile.Profile()
    t0 = time.perf_counter()
    pr.enable()
    weights, chi2_tot, chi2_kin, chi2_kinmap = solver.solve(
        orblib, ignore_existing_weights=True)
    pr.disable()
    print(f'\nsolve() {time.perf_counter()-t0:.1f}s, peak RSS {rss_gib():.1f} GiB')
    print(f'chi2_tot={chi2_tot:.6f} chi2_kin={chi2_kin:.6f} '
          f'sum(w)={np.nansum(weights):.10f} nonzero={int((weights>0).sum())}')
    pstats.Stats(pr).sort_stats('tottime').print_stats(20)
    np.save(f'/tmp/weights_{which}.npy', weights)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'adelie')
