import numpy as np
from . import dynamics

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
                 mod_dir=None,
                 settings=None):
        self.mod_dir = mod_dir
        self.settings = settings

    def solve(self):
        # self.write_executable()
        # self.execute()
        chi2, kinchi2 = self.read_output()
        return chi2, kinchi2

    def write_executable(self):
        # code to write the fortran executable here
        pass

    def execute(self):
        # code to run the executable here
        pass

    def read_output(self):
        ''' taken useful parts from triax_extract_chi2_iter in schw_domoditer,
        in particular lines 181-212
        '''
        # read amount of observables and kinematic moments
        fname = self.mod_dir + 'nn_kinem.out'
        a = self.__read_file_element(fname, [1, 1], [1, 2])
        ngh = np.int64(a[1])  # number of 'observables'
        nobs = np.int64(a[1])
        nvel = np.int64(a[0])
        ncon = np.int64(a[0])
        rows = 3 + np.arange(nobs)  # rows 1- 9
        cols = 3 + np.zeros(nobs, dtype=int)  # skip over text
        fname = self.mod_dir + 'nn_nnls.out'
        chi2vec = self.__read_file_element(fname, rows, cols)
        chi2vec = np.double(chi2vec)
        chi2 = sum(chi2vec)
        fname = self.mod_dir + 'nn_kinem.out'
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
