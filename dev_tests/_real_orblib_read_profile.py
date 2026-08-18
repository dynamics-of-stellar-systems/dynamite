"""Profile orblib.read_vel_histograms() on a real library.

It is the single most expensive step in a model once matrix assembly is fixed.

Run from ~/research/omegacen/dynamite_dataprep with the `main` conda python.
"""
import cProfile
import pstats
import time

import dynamite as dyn

CONFIG = 'NGC5139_config_veldist_combined_bigger.yaml'


def main():
    c = dyn.config_reader.Configuration(CONFIG, reset_logging=True)
    model = c.all_models.get_model_from_row(0)
    orblib = model.get_orblib()

    pr = cProfile.Profile()
    t0 = time.perf_counter()
    pr.enable()
    orblib.read_vel_histograms()
    pr.disable()
    print(f'\nread_vel_histograms: {time.perf_counter()-t0:.1f}s', flush=True)
    st = pstats.Stats(pr)
    print('\n=== by tottime ===')
    st.sort_stats('tottime').print_stats(22)


if __name__ == '__main__':
    main()
