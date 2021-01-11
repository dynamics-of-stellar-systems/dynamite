import subprocess
import myrand

seed = -4242

print("Comparison of portable random number generator myrand.py and"
      f" ran1_nr.f for seed={seed}")

print("Python random numbers from myrand.py:")
mr = myrand.MyRand(seed)
for i in range(10):
    print(mr.ran1())

print("Fortran random numbers from ran1_nr.f:")
subprocess.call(f"sed -i -e 's/^idum.*/idum = {seed}/' rantest.f90", shell=True)
compilestring = "gfortran -o rantest rantest.f90 ../legacy_fortran/ran1_nr.f"
runstring = "./rantest"
subprocess.call(f"{compilestring}\n{runstring}", shell=True)
print("Fortran version:")
subprocess.call("gfortran --version", shell=True)