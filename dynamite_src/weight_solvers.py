import os
import numpy as np
import subprocess
import dynamics

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
                 ml=None):
        self.system = system
        self.mod_dir = mod_dir
        self.settings = settings
        self.legacy_directory = legacy_directory
        self.ml=ml





    def solve(self):

        #set the current directory to the directory in which the models are computed
        cur_dir=os.getcwd()
        os.chdir(self.mod_dir)


        print("Fit the orbit library to the kinematic data.")

        cmdstr=self.write_executable()
        self.execute(cmdstr)

        #set the current directory to the dynamite directory
        os.chdir(cur_dir)

        chi2, kinchi2 = self.read_output()


        print("NNLS is finished.")

        return chi2, kinchi2

    def write_executable(self):

        '{:06.2f}'.format(self.ml)
        nn='ml'+'{:01.2f}'.format(self.ml)+'/nn'


        cmdstr = 'cmdd' + str(self.system.name) + '_' + str(int(np.random.uniform(0, 1) * 100000.0))

        txt_file = open(cmdstr, "w")
        txt_file.write('#!/bin/bash' + '\n')
        txt_file.write('# if the gzipped orbit library exist unzip it' + '\n')
        txt_file.write('test -e ../datfil/orblib.dat || bunzip2 -k  datfil/orblib.dat.bz2 ' + '\n')
        txt_file.write('test -e ../datfil/orblibbox.dat || bunzip2 -k  datfil/orblibbox.dat.bz2' + '\n')
        txt_file.write('test -e ' + str(nn) + '_kinem.out || ' +
                           self.legacy_directory +'/triaxnnls_CRcut < ' + str(nn) + '.in >>' +str(nn) + 'ls.log' + '\n') #TODO: specify which nnls to use
        txt_file.write('rm datfil/orblib.dat' + '\n')
        txt_file.write('rm datfil/orblibbox.dat' + '\n')
        txt_file.close()

        return cmdstr


    def execute(self,cmdstr):

          p4 = subprocess.call('bash ' + cmdstr, shell=True)

    def read_output(self):
        ''' taken useful parts from triax_extract_chi2_iter in schw_domoditer,
        in particular lines 181-212
        '''
        # read amount of observables and kinematic moments
        fname = self.mod_dir + 'ml'+'{:01.2f}'.format(self.ml)+'/nn_kinem.out'
        a = self.__read_file_element(fname, [1, 1], [1, 2])
        ngh = np.int64(a[1])  # number of 'observables'
        nobs = np.int64(a[1])
        nvel = np.int64(a[0])
        ncon = np.int64(a[0])
        rows = 3 + np.arange(nobs)  # rows 1- 9
        cols = 3 + np.zeros(nobs, dtype=int)  # skip over text
        fname = self.mod_dir + 'ml'+'{:01.2f}'.format(self.ml)+ '/nn_nnls.out'
        chi2vec = self.__read_file_element(fname, rows, cols)
        chi2vec = np.double(chi2vec)
        chi2 = sum(chi2vec)
        fname = self.mod_dir + 'ml'+'{:01.2f}'.format(self.ml)+'/nn_kinem.out'  #why is that in here twice?
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
