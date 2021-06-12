#!/bin/bash

testdir="test_all"

if [ $# -eq 1 ] && ( [ $1 = "?" ] || [ $1 = "-h" ] )
then
cat <<EOI
### Run a grid of DYNAMITE test scripts.

### A set of DYNAMITE test scripts are run with a variation
### of weight solvers and parameter generators, in total
### #(test scripts) * #(weight solvers/nnls solvers)
### * #(parameter generators) scenarios.

### The results are written into directory $testdir. Inside
### $testdir, each scenario is written into an individual
### subdirectory. Existing subdirectories are overwritten.
### Console output is captured in the directories' output.txt.

### HOW TO USE:
### 1. Verify that the test scripts and their config files are
###    in test_scripts and config_files, respectively.
### 2. Verify that the desired weight solver settings are in
###    weight_solvers and nnls_solvers, respectively.
### 3. Verify that the desired parameter generators are in
###    parameter_generators.
### 4. Verify that DYNAMITE has been installed correctly
###    ('python setup.py install --user'). If using SLURM,
###    make sure the python version matches the module loaded
###    in script variable slurm1.
### 4. Execute with $0 (run on local machine)
###    or $0 SLURM (use SLURM).

### '$0 ?' or '$0 -h' display this help message.
EOI
exit
fi

## --------------------
## Test scripts: enter desired test scripts/config file combinations here.
## Note that config_files must contain matching configuration files.
test_scripts=('test_nnls.py' 'test_orbit_losvds.py')
config_files=('user_test_config_ml.yaml' 'user_test_config.yaml')

## --------------------
## Weight solvers: enter desired type/nnls_solver combinations here.
## Note that nnls_solvers must contain matching nnls_solver entries
## for the configuration files.
weight_solvers=('LegacyWeightSolver' 'NNLS' 'NNLS')
nnls_solvers=('1' 'scipy' 'cvxopt')

## --------------------
## Parameter generators: enter desired generator_types here
parameter_generators=('LegacyGridSearch' 'GridWalk' 'FullGrid')

## --------------------
## SLURM script: slurm0 and slurm1 are used to build SLURM batch scripts
## used for running each test as an individual batch job.
## ** Activated with command-line argument SLURM **
slurm0="#!/bin/bash \n#\n"
slurm0="${slurm0}#SBATCH --qos=p71474_0096\n"
slurm0="${slurm0}#SBATCH --job-name="
slurm1="\n#SBATCH -N 1\n"
slurm1="${slurm1}#SBATCH --output=\"dyn_%%j.out\"\n"
slurm1="${slurm1}#SBATCH --error=\"dyn_%%j.err\"\n"
slurm1="${slurm1}#SBATCH --ntasks-per-node=16\n"
slurm1="${slurm1}#SBATCH --ntasks-per-core=1\n"
slurm1="${slurm1}pwd; hostname; date\n"
slurm1="${slurm1}\nmodule purge\n"
slurm1="${slurm1}module load python/3.8.0-gcc-9.1.0-6jpq4wd\n"
slurm1="${slurm1}module load py-setuptools\n\n"
joblist=""

printf "Starting $0, $(date)\n"

PY=python

cwd=$(pwd)
[[ -d "$testdir" ]] || mkdir "$testdir"
total_count=$((${#weight_solvers[*]}*${#parameter_generators[*]}*${#test_scripts[*]}))
typeset -i ws_i=0 ws_n=${#weight_solvers[*]} script_i script_n=${#test_scripts[*]} count=1
while (( ws_i < ws_n ))
do
  ws=${weight_solvers[ws_i]}
  nnls=${nnls_solvers[ws_i]}
  for pg in ${parameter_generators[*]}
  do
    script_i=0
    while (( script_i < script_n ))
    do
      script=${test_scripts[script_i]}
      config=${config_files[script_i]}
      folder="${testdir}/${script_i}_${ws}_${nnls}_$pg"
      # create a fresh test folder
      rm -rf $folder
      mkdir $folder
      # copy sript and comparison data to test folder
      cp $script $folder
      cp -r data $folder
      # adapt and copy configuration file to folder
      awk 'BEGIN{ws=0} {if(match($1, "weight_solver_settings")){ws=1};if(ws && match($1,"type")){$0="    type: \"'$ws'\"";ws=0};print}' $config | \
      sed -e "s/    nnls_solver:.*/    nnls_solver: \"$nnls\"/" \
          -e "s/    generator_type:.*/    generator_type: \"$pg\"/" \
          -e 's/    input_directory: "\(.*\)"/    input_directory: "..\/..\/\1"/' > $folder/$config
      cd $folder
      if [ $# -eq 1 ] && [ $1 = "SLURM" ]
      then
        jobname=$(echo $folder | cut -d/ -f2)
        printf "${slurm0}${jobname}${slurm1}$PY $script &> output.txt\ndate" > run.slrm
        jobid=$(sbatch run.slrm | cut -d' ' -f4)
        printf "$count of $total_count: SLURM job executing $script using $ws/$nnls and ${pg}: $jobid\n"
        if [ -z "$joblist" ]
        then
          joblist=$jobid
        else
          joblist="${joblist},$jobid"
        fi
      else
        printf "$count of $total_count: Executing $script using $ws/$nnls and $pg... "
        # run the script
        $PY $script &> output.txt
        exitcode=$?
        if [ $exitcode -eq 0 ]
        then
          printf "OK"
        else
          printf "ERROR"
        fi
        printf ".\n"
      fi
      cd $cwd
      count+=1
      script_i+=1
    done
  done
  ws_i+=1
done

printf "$0 done, $(date)\n"
if [ $# -eq 1 ] && [ $1 = "SLURM" ]
then
  jobstat="sacct --jobs=$joblist --format=JobID%15,JobName%40,Start,Submit,NodeList,NCPUS,State,ExitCode"
  $jobstat
  printf "\nDon't worry if the table abve does not show all jobs. Check job status with:\n\n"
  echo ${jobstat}
  echo
fi
