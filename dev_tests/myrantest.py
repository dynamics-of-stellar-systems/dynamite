#!/usr/bin/env python3
import sys
import subprocess
import logging
import numpy as np
import dynamite.myrand as myrand
from dynamite.config_reader import DynamiteLogging

SEED = -4242

def run_random_number_test(ran_seed=SEED, n_ran=10, make_comp=False):
    logger.info("Comparison of portable random number generator myrand.py and"
                f" ran1_nr.f for seed={ran_seed}")

    rtol = 0.
    atol = 1e-16

    fname = f"data/randata{ran_seed}.txt"
    if not make_comp:
        ran_saved = np.loadtxt(fname)
        n_ran = len(ran_saved)
        logger.debug(f"{n_ran} random numbers read")

    mr = myrand.MyRand(ran_seed)
    ran_P = [mr.ran1() for i in range(n_ran)]
    logger.info(f"Generated {n_ran} Python random numbers from myrand.py")
    logger.debug(ran_P)

    subprocess.run(f"sed -i -e 's/^idum.*/idum = {ran_seed}/' rantest.f90", \
                   shell=True, check=True)
    subprocess.run(f"sed -i -e 's/^do i.*/do i=1,{n_ran-1}/' rantest.f90", \
                   shell=True, check=True)
    compilestring="gfortran -o rantest rantest.f90 ../legacy_fortran/ran1_nr.f"
    runstring = "./rantest"
    out_F = subprocess.run(f"{compilestring}\n{runstring}", \
                           capture_output=True, shell=True, check=True)
    ran_F = [float(n) for n in out_F.stdout.split()]
    logger.info(f"Generated {n_ran} Fortran random numbers from ran1_nr.f")
    logger.debug(ran_F)

    is_same = np.allclose(ran_P, ran_F, rtol=rtol, atol=atol)
    if is_same:
        logger.info(f"All {n_ran} generated random number pairs match")
    else:
        text = f"Mismatch of at leat one of the {n_ran} generated random " \
               "number pairs!"
        logger.error(text)
        ValueError(text)

    if not make_comp:
        is_same = np.allclose(ran_P, ran_saved, rtol=rtol, atol=atol)
        if is_same:
            logger.info(f"SUCCESS - all {n_ran} random numbers match the "
                        f"saved values (file: {fname})")
        else:
            logger.info(f"Mismatch of at leat one of the {n_ran} random "
                        f"numbers (file: {fname})!")
    else:
        np.savetxt(fname, ran_P)
        logger.info(f"{n_ran} random numbers saved in file {fname}")

    version_P = sys.version.split(sep='\n')[0]
    logger.debug(f"Python version: {version_P}")
    info_F = subprocess.run("gfortran --version", capture_output=True, \
                            shell=True, check=True)
    version_F = info_F.stdout.decode('utf-8').split(sep='\n')[0]
    logger.debug(f"Fortran version: {version_F}")

if __name__ == '__main__':
    logger = logging.getLogger()
    DynamiteLogging(logfile='myrantest.log')
    run_random_number_test()

# end
