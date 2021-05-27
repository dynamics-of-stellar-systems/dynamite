#!/bin/bash

### Run a grid of test scripts

### This test script executes the scripts given in test_scripts.
### In the configuration files, the weight solvers and parameter
### generators are varied according to the entries in
### weight_solvers/nnls_solvers and parameter_generators,
### respectively.

### The results are written into $testdir. Inside this directory, a new
### directory holding all the test data is created for each combination
### of <test script> - <weight solver> - <parameter generator>.
### The console output is captured in the directories' output.txt.

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

printf "Starting $0, $(date)\n"

PY=python

cwd=$(pwd)
testdir="test_all"
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
      folder="${testdir}/test_${script_i}_${ws}_${nnls}_$pg"
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
      printf "$count of $total_count: Executing $script using $ws/$nnls and $pg... "
      cd $folder
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
      count+=1
      script_i+=1
    done
  done
  ws_i+=1
done

printf "$0 done, $(date)\n"
