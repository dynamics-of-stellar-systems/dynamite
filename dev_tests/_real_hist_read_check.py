"""Fingerprint the velocity histograms read from a real orbit library.

Run once on each git revision and compare: the bulk parsers must reproduce the
per-record FortranFile loop exactly. Uses sha256 of the raw bytes rather than
holding two copies of a ~10 GiB array in memory.

Run from ~/research/omegacen/dynamite_dataprep with the `main` conda python:
    python .../_real_hist_read_check.py
"""
import hashlib
import time

import dynamite as dyn

CONFIG = 'NGC5139_config_veldist_combined_bigger.yaml'


def main():
    c = dyn.config_reader.Configuration(CONFIG, reset_logging=True)
    model = c.all_models.get_model_from_row(0)
    orblib = model.get_orblib()

    t0 = time.perf_counter()
    orblib.read_vel_histograms()
    dt = time.perf_counter() - t0
    print(f'\nread_vel_histograms: {dt:.1f}s', flush=True)

    for i, h in enumerate(orblib.vel_histograms):
        y = h.y
        digest = hashlib.sha256(y.tobytes()).hexdigest()[:32]
        print(f'  hist {i}: shape={y.shape} dtype={y.dtype} '
              f'nonzero={int((y != 0).sum())}')
        print(f'           sum={y.sum()!r}')
        print(f'           sha256[:32]={digest}')


if __name__ == '__main__':
    main()
