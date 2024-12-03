#!/bin/bash

testdir="../docs/tutorial_notebooks/"

if [ $# -eq 1 ] && ( [ $1 = "?" ] || [ $1 = "-h" ] )
then
cat <<EOI
### Run the DYNAMITE tutorial notebooks from the command line

### HOW TO USE:
### 1. Check that you are in dev_tests, or adjust the 'testdir'
### variable to point to a different location for the notebooks
### 2. Run ./test_notebooks.sh
### '$0 ?' or '$0 -h' display this help message.
EOI
exit
fi

notebooks="1_data_prep_for_gauss_hermite.ipynb 2_quickstart.ipynb 3_model_iterations_and_plots.ipynb 4_BayesLOSVD.ipynb 5_parameter_space.ipynb 6_orbits_and_weights.ipynb 7_orbital_distributions.ipynb"

cd $testdir

for n in $notebooks; do
    jupyter execute $n >/dev/null 2>&1 && echo Passed on notebook: $n || echo Failed on notebook: $n
done
