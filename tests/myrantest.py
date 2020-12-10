import subprocess
import myrand

seed = -4242

print(f"Python random numbers for seed = {seed}:")
mr = myrand.MyRand(seed)
for i in range(10):
    print(mr.ran1())

print(f"Fortran random numbers for seed = {seed}:")
subprocess.call(f"sed -i -e 's/^idum.*/idum = {seed}/' rantest.f90", shell=True)
compilestring = "gfortran -o rantest rantest.f90 ../legacy_fortran/ran1_nr.f"
runstring = "./rantest"
subprocess.call(f"{compilestring}\n{runstring}", shell=True)
