import os
import numpy as np
import subprocess
import orblib
from astropy import table

class WeightSolver(object):

    def __init__(self,
                 weight_solver_args={}
                 ):
        self.weight_solver_args = weight_solver_args

    def set_kinematics(self, system):
        self.kinematics = []
        for component in system.cmp_list:
            self.kinematics += component.kinematic_data

    def set_orb_lib(self, orb_lib):
        self.orb_lib = orb_lib

    def get_observables_from_orbits(self):
        self.orb_observed = []
        for kin_data0 in self.kinematics:
            orb_obs0 = kin_data0.transform_orbits_to_observables(self.orb_lib)
            self.orb_observed += [orb_obs0]

    def solve(self):
        # placeholder function to solve for weights given
        # self.kinematics.values and self.orb_observed
        # return wts, chi2
        return 0, 0


class NNLS(WeightSolver):

    def __init__(self,
                 **kwargs):
        super(NNLS, self).__init__(*kwargs)     # initialise parent class

    def solve(self):
        # actual code to do NNLS
        # return MAP_weight
        return 0, 0


class LegacyWeightSolver(WeightSolver):

    def __init__(self,
                 system=None,
                 mod_dir=None,
                 settings=None,
                 legacy_directory=None,
                 ml=None,
                 executor=None):
        self.system = system
        self.mod_dir = mod_dir
        self.settings = settings
        self.legacy_directory = legacy_directory
        self.ml=ml
        self.executor = executor
        self.mod_dir_with_ml = self.mod_dir + 'ml' + '{:01.2f}'.format(self.ml)
        self.fname_nn_kinem = self.mod_dir_with_ml + '/nn_kinem.out'
        self.fname_nn_nnls = self.mod_dir_with_ml + '/nn_nnls.out'

    def solve(self):
        check1 = os.path.isfile(self.fname_nn_kinem)
        check2 = os.path.isfile(self.fname_nn_nnls)
        if not check1 or not check2:
            # set the current directory to the directory in which the models are computed
            cur_dir = os.getcwd()
            os.chdir(self.mod_dir)
            print("Fit the orbit library to the kinematic data.")
            cmdstr = self.executor.write_executable_for_weight_solver(self.ml)
            self.executor.execute(cmdstr)
            #set the current directory to the dynamite directory
            os.chdir(cur_dir)
            print("NNLS problem solved")
        else:
            print("NNLS solution read from existing output")
        weights = self.read_weights()
        chi2, kinchi2 = self.read_chi2()
        return chi2, kinchi2

    def read_weights(self):
        fname = self.mod_dir_with_ml + '/nn_orb.out'
        col_names = ['orb_idx',
                     'E_idx',
                     'I2_idx',
                     'I3_idx',
                     'I_dont_know', # NOTE: this column = 0 - unclear what it is
                     'orb_type',
                     'weight',
                     'lcut'] # NOTE: see lines 1321-1322 of triaxnnls_CRcut.f90
        # NOTE: possible that column 'lcut' is not present when a different
        # "triaxnnls" file is used
        dtype = [int, int, int, int, int, int, np.float64, int]
        weights = np.genfromtxt(fname,
                                skip_header=1,
                                names=col_names,
                                dtype=dtype)
        weights = table.Table(weights)
        # TODO: find out what the following column should be doing
        weights.remove_column('I_dont_know')
        self.weights = weights

    def read_chi2(self):
        ''' taken useful parts from triax_extract_chi2_iter in schw_domoditer,
        in particular lines 181-212
        '''
        # read amount of observables and kinematic moments
        fname = self.fname_nn_kinem
        a = self.__read_file_element(fname, [1, 1], [1, 2])
        ngh = np.int64(a[1])  # number of 'observables'
        nobs = np.int64(a[1])
        nvel = np.int64(a[0])
        ncon = np.int64(a[0])
        rows = 3 + np.arange(nobs)  # rows 1- 9
        cols = 3 + np.zeros(nobs, dtype=int)  # skip over text
        fname = self.fname_nn_nnls
        chi2vec = self.__read_file_element(fname, rows, cols)
        chi2vec = np.double(chi2vec)
        chi2 = sum(chi2vec)
        fname = self.fname_nn_kinem
        ka = np.genfromtxt(fname, skip_header=1)
        k = np.arange(ngh) * 3 + 3
        kinchi2 = sum(sum(pow(((ka[:, k] - ka[:, k + 1]) / ka[:, k + 2]), 2.0)))
        return chi2, kinchi2

    def __read_file_element(self, infile, rows, cols):
        """Taken from schw_misc
        !@brief read fields in a tabular data according to the their row/column.
        @details Function description.
        @param[in] infile input file
        @param[in] row array of locations (row), indexing starts from 1.
        @param[in] cols array of locations (column), indexing starts from 1.
        """
        lines = [line.rstrip('\n').split() for line in open(infile)]
        output=[]
        for i in range(0, len(rows)):
            output.append(lines[rows[i] - 1][cols[i] - 1])
        return output

# end
