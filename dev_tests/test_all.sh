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

### Inspect the code and its annotations for further details.

### '$0 ?' or '$0 -h' display this help message.
EOI
exit
fi

## --------------------
## Test scripts: enter desired test scripts/config file combinations here.
## Note that config_files must contain element-wise matching configuration
## files and ncores the number of CPU to be used for each test script.
test_scripts=('test_nnls.py')
config_files=('user_test_config_ml.yaml')
ncores=('8')

## --------------------
## Weight solvers: enter desired type/nnls_solver combinations here.
## Note that nnls_solvers must contain element-wise matching nnls_solver
## entries for the respective configuration files.
weight_solvers=('LegacyWeightSolver' 'NNLS')
nnls_solvers=('1' 'scipy')

## --------------------
## Parameter generators: enter desired generator_types here
parameter_generators=('LegacyGridSearch' 'GridWalk' 'FullGrid')

## --------------------
## If using a cluster and Slurm, the variable slurm_array determines
## which 'dimension' will be used for the Slurm job array:
## "script" will use as many nodes as there are test_scripts
## "ws" will use as many nodes as weight_solvers
## "pg" will use as many nodes as parameter_generators
slurm_array="script"  # create one array job per "script", "ws", or "pg"

## --------------------
## SLURM script: slurm0 and slurm1 are used to build SLURM batch scripts
## used for running each test as an individual batch job.
## ** Activated with command-line argument SLURM **
slurm0="#!/bin/bash \n#\n"
slurm0="${slurm0}#SBATCH --qos=p71474_0096\n"
slurm0="${slurm0}#SBATCH --job-name=DYNAMITE_test\n"
slurm0="${slurm0}#SBATCH -N 1\n"
slurm0="${slurm0}#SBATCH --output=\"dyn_%%j.out\"\n"
slurm0="${slurm0}#SBATCH --error=\"dyn_%%j.err\"\n"
slurm0="${slurm0}#SBATCH --ntasks-per-node=96\n"
slurm0="${slurm0}#SBATCH --ntasks-per-core=2\n"
slurm0="${slurm0}#SBATCH --array=0-"
slurm1="\n\npwd; hostname; date\n"
slurm1="${slurm1}\nmodule purge\n"
slurm1="${slurm1}module load python/3.8.0-gcc-9.1.0-6jpq4wd\n"
slurm1="${slurm1}module load py-setuptools\n\n"

### Script start

printf "Starting $0, $(date)\n"

PY=python

cwd=$(pwd)
[[ -d "$testdir" ]] || mkdir "$testdir"
total_count=$((${#weight_solvers[*]}*${#parameter_generators[*]}*${#test_scripts[*]}))
typeset -i ws_i=0 ws_n=${#weight_solvers[*]} script_i script_n=${#test_scripts[*]} count=1
typeset -i pg_i pg_n=${#parameter_generators[*]}

# Create Slurm script header
if [ $# -eq 1 ] && [ $1 = "SLURM" ]
then
  max_array_var=${slurm_array}_n
  i_array_var=${slurm_array}_i
  run=run$(date +%s).slrm
  printf "$slurm0$((${!max_array_var}-1))$slurm1" > $testdir/$run
fi

# Main loop
while (( ws_i < ws_n ))
do
  ws=${weight_solvers[ws_i]}
  nnls=${nnls_solvers[ws_i]}
  pg_i=0
  for pg in ${parameter_generators[*]}
  do
    script_i=0
    while (( script_i < script_n ))
    do
      script=${test_scripts[script_i]}
      config=${config_files[script_i]}
      folder="${script_i}_${ws}_${nnls}_$pg"
      # create a fresh test folder
      rm -rf $testdir/$folder
      mkdir $testdir/$folder
      # copy script and comparison data to test folder
      cp $script $testdir/$folder
      cp -r data $testdir/$folder
      # adapt and copy DYNAMITE configuration file to test folder
      # first, set fixed: False for the black hole mass (must be the first
      # parameter in the bh component with existing fixed: ... entry)
      awk 'BEGIN{m=0} {if(match($1, "bh:")){m=1};if(m && match($1,"fixed:")){$0="                fixed: False";m=0};print}' $config | \
      # then, inject the weight solver into the configuration file
      awk 'BEGIN{ws=0} {if(match($1, "weight_solver_settings")){ws=1};if(ws && match($1,"type")){$0="    type: \"'$ws'\"";ws=0};print}' | \
      # at last, inject the nnls solver, generator type, input directory, and ncpus
      sed -e "s/    nnls_solver:.*/    nnls_solver: \"$nnls\"/" \
          -e "s/    generator_type:.*/    generator_type: \"$pg\"/" \
          -e 's/    input_directory: "\(.*\)"/    input_directory: "..\/..\/\1"/' \
          -e "s/    ncpus:.*/    ncpus: ${ncores[script_i]}/" > $testdir/$folder/$config
      if [ $# -eq 1 ] && [ $1 = "SLURM" ]
      then
        # Add scenario to Slurm file and sort into job array
        printf "if [ \$SLURM_ARRAY_TASK_ID -eq ${!i_array_var} ]\nthen\n  cd $folder && $PY $script &> output.txt &\nfi\n" >> $testdir/$run
      else
        # Execute local jobs
        cd $testdir/$folder
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
        cd $cwd
      fi
      count+=1
      script_i+=1
    done
  pg_i+=1
  done
  ws_i+=1
done

printf "$0 done, $(date)\nLook into $testdir.\n"

if [ $# -eq 1 ] && [ $1 = "SLURM" ]
then
  cd $testdir
  # Finish writing Slurm script
  printf "\nwait\n" >> $run
  # Submit Slurm jobs
  jobid=$(sbatch $run | cut -d' ' -f4)
  cd ..
  printf "Submitted SLURM job $jobid comprising $total_count scenarios in array of ${!max_array_var}:\n"
  jobstat="sacct --jobs=$jobid --format=JobID%15,JobName%15,Submit,Start,End,NodeList,State,ExitCode"
  squeue="squeue --user `whoami`"
  echo
  $squeue
  echo
  $jobstat
  echo
  printf "\nDon't worry if the table abve doesn't show the job yet. Check job status with:\n\n"
  echo ${squeue}
  echo or
  echo ${jobstat}
fi
