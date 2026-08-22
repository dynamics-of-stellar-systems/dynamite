#!/usr/bin/env python3
"""Check that integrating an orbit library in chunks changes nothing.

Builds the same orbit library several times through DYNAMITE, varying only
``multiprocessing_settings: orblib_chunks``, and requires that the merged
libraries are byte-identical to the single-process one and that the weight
solver returns identical chi2.

Each configuration runs in its own copy of the model tree, so they can run
concurrently without contending over ``datfil``. On APFS the copies are clones
and cost nothing.

Chunk count 7 is included deliberately: it does not divide the orbit count, so
it exercises the distribution of the remainder.

Run it at dithering > 1 as well as dithering 1. The two differ in how many
starting points the Fortran integrates per orbit bundle, and a bug in that
count is invisible at dithering 1, where the correction factor is 1:

    python test_orblib_chunking.py --config user_test_config_ml_dither2.yaml

That needs a model tree for the same configuration to exist first (this test
reuses begin.dat rather than spending minutes regenerating it), so seed one
with a single-model run of the same config before the first use. Do not point
it at a tree built from a *different* config: the grid in begin.dat would not
match and the Fortran would refuse the chunk ranges.

Usage:
    python test_orblib_chunking.py [--config user_test_config_ml.yaml]
"""

import argparse
import filecmp
import glob
import os
import shutil
import subprocess
import sys
import time
from concurrent import futures

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import dynamite as dyn                                    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
LIB_FILES = ('orblib_qgrid.dat.bz2', 'orblib_losvd_hist.dat.bz2',
             'orblibbox_qgrid.dat.bz2', 'orblibbox_losvd_hist.dat.bz2',
             # analysis.py and plotter.py read these for lambda_z, and
             # model.py counts them towards a complete library
             'orblib.dat_orbclass.out', 'orblibbox.dat_orbclass.out')


def make_work_dir(root, name, config, out_dir, in_dir, n_chunks):
    """Copy the model tree and write a config with the given chunk count."""
    work = os.path.join(root, name)
    if os.path.exists(work):
        shutil.rmtree(work)
    os.makedirs(work)
    for d in (out_dir, in_dir):
        src = os.path.join(HERE, d)
        if os.path.isdir(src):
            # -c asks for APFS clones; harmless where unsupported
            subprocess.run(['cp', '-Rc', src, os.path.join(work, d)],
                           check=True)
    text = open(os.path.join(HERE, config)).read()
    if 'orblib_chunks' in text:
        raise ValueError(f'{config} already sets orblib_chunks; this test '
                         'needs to control it')
    text = text.replace('multiprocessing_settings:',
                        'multiprocessing_settings:\n'
                        f'    orblib_chunks: {n_chunks}', 1)
    cfg = os.path.join(work, 'cfg.yaml')
    open(cfg, 'w').write(text)
    return work, cfg


def build_subprocess(work, cfg, out_dir):
    """Run one build in its own process.

    Each build changes directory and reconfigures DYNAMITE's logging, both of
    which are process-global, so builds cannot share an interpreter.

    Returns
    -------
    tuple
        ``(datfil, seconds, chi2, kinchi2)``

    """
    r = subprocess.run(
        [sys.executable, os.path.abspath(__file__), '--build', work, cfg,
         out_dir],
        capture_output=True, text=True)
    for line in reversed(r.stdout.splitlines()):
        if line.startswith('RESULT '):
            datfil, secs, chi2, kinchi2 = line[7:].rsplit('\t', 3)
            return datfil, float(secs), float(chi2), float(kinchi2)
    raise RuntimeError(f'build in {work} produced no result:\n'
                       f'{r.stdout[-2000:]}\n{r.stderr[-2000:]}')


def build(work, cfg, out_dir):
    """Force a fresh integration in ``work``; return files, time and chi2.

    Works from an empty output directory: the model tree and the initial
    conditions are built if absent. Reusing a tree left over from a *different*
    configuration is how this test came to run against a ``begin.dat`` for the
    wrong grid, so it no longer requires one to exist.
    """
    cwd = os.getcwd()
    os.chdir(work)
    try:
        t_0 = time.perf_counter()
        c = dyn.config_reader.Configuration(os.path.basename(cfg),
                                            reset_logging=True,
                                            reset_existing_output=False)
        if len(c.all_models.table) > 0:
            mod = c.all_models.get_model_from_row(0)
        else:
            from astropy import table as aptable
            vals = {p.name: [p.par_value] for p in c.parspace}
            parset = aptable.Table(
                {n: vals[n] for n in c.parspace.par_names})[0]
            mod = dyn.model.Model(config=c, parset=parset)
        # discard any existing library so the integration really runs, but keep
        # begin.dat/beginbox.dat: regenerating the initial conditions is slow
        # and they do not depend on the chunk count
        datfil = os.path.join(work, mod.directory_noml, 'datfil')
        # normally created when the model iterator assigns directories
        os.makedirs(datfil, exist_ok=True)
        os.makedirs(os.path.join(work, mod.directory_noml, 'infil'),
                    exist_ok=True)
        for f in glob.glob(os.path.join(datfil, 'orblib*')):
            os.remove(f)
        for f in ('tube_done', 'box_done', 'tube_box_done', 'building.lock'):
            if os.path.isfile(os.path.join(datfil, f)):
                os.remove(os.path.join(datfil, f))
        orblib = mod.get_orblib()
        elapsed = time.perf_counter() - t_0
        ws = dyn.weight_solvers.NNLS(config=c, model=mod)
        _, chi2, kinchi2, _ = ws.solve(orblib, ignore_existing_weights=True)
    finally:
        os.chdir(cwd)
    return datfil, elapsed, float(chi2), float(kinchi2)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == '--build':
        # internal single-build mode, invoked by build_subprocess
        datfil, secs, chi2, kinchi2 = build(*sys.argv[2:5])
        print(f'RESULT {datfil}\t{secs}\t{chi2}\t{kinchi2}')
        return 0

    ap = argparse.ArgumentParser()
    ap.add_argument('--config', default='user_test_config_ml.yaml')
    ap.add_argument('--chunks', default='1,2,4,7,8')
    ap.add_argument('--workdir', default='/tmp/dyn_chunk_test')
    args = ap.parse_args()

    counts = [int(x) for x in args.chunks.split(',')]
    if 1 not in counts:
        counts.insert(0, 1)

    cfg_text = open(os.path.join(HERE, args.config)).read()

    def setting(name):
        for line in cfg_text.splitlines():
            if line.strip().startswith(name):
                return line.split(':', 1)[1].strip().strip('"\'').rstrip('/')
        raise KeyError(name)

    out_dir, in_dir = setting('output_directory'), setting('input_directory')

    print(f'dynamite {dyn.__version__} from {dyn.__path__[0]}')
    print(f'config {args.config}, chunk counts {counts}', flush=True)

    os.makedirs(args.workdir, exist_ok=True)
    jobs = {n: make_work_dir(args.workdir, f'chunks{n}', args.config,
                             out_dir, in_dir, n) for n in counts}

    results = {}
    with futures.ThreadPoolExecutor(len(jobs)) as ex:
        fut = {ex.submit(build_subprocess, w, c, out_dir): n
               for n, (w, c) in jobs.items()}
        for f in futures.as_completed(fut):
            n = fut[f]
            results[n] = f.result()
            print(f'  chunks={n}: {results[n][1]:6.1f}s  '
                  f'chi2={results[n][2]:.2f}  kinchi2={results[n][3]:.2f}',
                  flush=True)

    ref_datfil, ref_t, ref_chi2, ref_kin = results[1]
    print(f'\n{"chunks":>7s} {"files match":>12s} {"chi2 match":>11s} '
          f'{"time":>9s}')
    ok = True
    for n in counts:
        datfil, elapsed, chi2, kinchi2 = results[n]
        same_files = all(
            filecmp.cmp(os.path.join(ref_datfil, f), os.path.join(datfil, f),
                        shallow=False) for f in LIB_FILES)
        same_chi2 = (chi2 == ref_chi2) and (kinchi2 == ref_kin)
        ok &= same_files and same_chi2
        print(f'{n:7d} {str(same_files):>12s} {str(same_chi2):>11s} '
              f'{elapsed:8.1f}s')

    print(f'\nchunked == unchunked: {"PASS" if ok else "FAIL"}')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
