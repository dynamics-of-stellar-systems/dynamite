#!/bin/bash

# Executes all .py files in current directory
# Intended for running all DYNAMITE dev_tests

PY=python

allgood=TRUE

for script in *.py
do
  name=${script%.*}
  outfile="output_${name}.txt"
  printf "Running $script ... "
  starttime=$(date +%s)
  # Execute test script and redirect all output into outfile
  $PY $script &> $outfile
  exitcode=$?
  endtime=$(date +%s)
  duration=$(bc <<< "$endtime - $starttime")
  if [ $exitcode -eq 0 ]
  then
    printf "OK"
  else
    allgood=FALSE
    printf "ERROR"
  fi
  printf ". Runtime ${duration}s. Look at ${outfile}"
  # delete previously saved logs
  rm -f *_${name}.log
  # find current logs and preserve them by renaming
  for logfile in $(ls *.log)
  do
    modtime=$(date -r $logfile +%s)
    if [ $modtime -gt $starttime ]
    then
      newlog="${logfile%.*}_${name}.log"
      mv $logfile $newlog
      printf ", $newlog"
    fi
  done
  printf ".\n" 
done

echo "***** Files to look at (raw output from test scripts):"
grep Look output_*.txt

[ $allgood = "TRUE" ] && echo "**SUCCESS** - all tests executed"
[ $allgood = "FALSE" ] && echo "At least one test **FAILED** - check .txt and .log files"
