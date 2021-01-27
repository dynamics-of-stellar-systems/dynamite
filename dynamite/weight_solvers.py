import os
import numpy as np
from astropy import table
import subprocess
from scipy import optimize

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
        self.mod_dir_with_ml = self.mod_dir + 'ml' + '{:01.2f}'.format(self.ml)
        self.fname_nn_kinem = self.mod_dir_with_ml + '/nn_kinem.out'
        self.fname_nn_nnls = self.mod_dir_with_ml + '/nn_nnls.out'
        # prepare fortran input file for nnls
        self.create_fortran_input_nnls(self.mod_dir, ml)

    def create_fortran_input_nnls(self,path,ml):

        #for the ml the model is only scaled. We therefore need to know what is the ml that was used for the orbit library
        infile=path+'infil/parameters.in'
        lines = [line.rstrip('\n').split() for line in open(infile)]
        ml_orblib=float((lines[-9])[0])

        #-------------------
        #write nn.in
        #-------------------

        text='infil/parameters.in' +'\n' + \
        str(self.settings['regularisation'])   + '                                  [ regularization strength, 0 = no regularization ]' +'\n'  + \
        'ml'+ '{:01.2f}'.format(ml) + '/nn' +'\n' + \
        'datfil/mass_qgrid.dat' +'\n' + \
        'datfil/mass_aper.dat' +'\n' + \
        str(self.settings['number_GH']) + '	                           [ # of GH moments to constrain the model]' +'\n' + \
        'infil/kin_data.dat' +'\n' + \
        str(self.settings['GH_sys_err']) + '    [ systemic error of v, sigma, h3, h4... ]' + '\n' + \
        str(self.settings['lum_intr_rel_err']) + '                               [ relative error for intrinsic luminosity ]' +'\n' + \
        str(self.settings['sb_proj_rel_err']) + '                               [ relative error for projected SB ]' + '\n' + \
        str(np.sqrt(ml/ml_orblib))  + '                                [ scale factor related to M/L, sqrt( (M/L)_k / (M/L)_ref ) ]' + '\n' + \
        f'datfil/orblib_{ml}.dat' +'\n' + \
        f'datfil/orblibbox_{ml}.dat' +'\n' + \
        str(self.settings['nnls_solver']) + '                                  [ nnls solver ]'

        nn_file= open(path+'ml'+'{:01.2f}'.format(ml)+'/nn.in',"w")
        nn_file.write(text)
        nn_file.close()

    def solve(self, orblib, ret_weights=False):
        # note: orblib is provided as an argument for consistency with future
        # WeightSolver implementations but are not used in this legacy method
        check1 = os.path.isfile(self.fname_nn_kinem)
        check2 = os.path.isfile(self.fname_nn_nnls)
        if not check1 or not check2:
            # set the current directory to the directory in which the models are computed
            cur_dir = os.getcwd()
            os.chdir(self.mod_dir)
            print("Fit the orbit library to the kinematic data.")
            cmdstr = self.write_executable_for_weight_solver(self.ml)
            p = subprocess.call('bash '+cmdstr, shell=True)
            #set the current directory to the dynamite directory
            os.chdir(cur_dir)
            print("NNLS problem solved")
        else:
            print("NNLS solution read from existing output")
        weights = self.read_weights()
        chi2, kinchi2 = self.read_chi2()
        return chi2, kinchi2

    def write_executable_for_weight_solver(self, ml):
        nn = f'ml{ml:01.2f}/nn'
        cmdstr = f'cmd_nnls_{ml}'
        txt_file = open(cmdstr, "w")
        txt_file.write('#!/bin/bash' + '\n')
        txt_file.write('# if the gzipped orbit library exist unzip it' + '\n')
        txt_file.write(f'test -e datfil/orblib_{ml}.dat || bunzip2 -c  datfil/orblib.dat.bz2 > datfil/orblib_{ml}.dat' + '\n')
        txt_file.write(f'test -e datfil/orblibbox_{ml}.dat || bunzip2 -c  datfil/orblibbox.dat.bz2 > datfil/orblibbox_{ml}.dat' + '\n')
        txt_file.write('test -e ' + str(nn) + '_kinem.out || ' +
                       self.legacy_directory +
                       f'/triaxnnls_CRcut < {nn}.in >> {nn}ls.log' + '\n')
        #TODO: specify which nnls to use
        txt_file.write(f'rm datfil/orblib_{ml}.dat' + '\n')
        txt_file.write(f'rm datfil/orblibbox_{ml}.dat' + '\n')
        txt_file.close()
        return cmdstr

    def create_fortran_input_nnls(self,path,ml):

        #for the ml the model is only scaled. We therefore need to know what is the ml that was used for the orbit library
        infile=path+'infil/parameters.in'
        lines = [line.rstrip('\n').split() for line in open(infile)]
        ml_orblib=float((lines[-9])[0])

        #-------------------
        #write nn.in
        #-------------------

        text='infil/parameters.in' +'\n' + \
        str(self.settings['regularisation'])   + '                                  [ regularization strength, 0 = no regularization ]' +'\n'  + \
        'ml'+ '{:01.2f}'.format(ml) + '/nn' +'\n' + \
        'datfil/mass_qgrid.dat' +'\n' + \
        'datfil/mass_aper.dat' +'\n' + \
        str(self.settings['number_GH']) + '	                           [ # of GH moments to constrain the model]' +'\n' + \
        'infil/kin_data.dat' +'\n' + \
        str(self.settings['GH_sys_err']) + '    [ systemic error of v, sigma, h3, h4... ]' + '\n' + \
        str(self.settings['lum_intr_rel_err']) + '                               [ relative error for intrinsic luminosity ]' +'\n' + \
        str(self.settings['sb_proj_rel_err']) + '                               [ relative error for projected SB ]' + '\n' + \
        str(np.sqrt(ml/ml_orblib))  + '                                [ scale factor related to M/L, sqrt( (M/L)_k / (M/L)_ref ) ]' + '\n' + \
        f'datfil/orblib_{ml}.dat' +'\n' + \
        f'datfil/orblibbox_{ml}.dat' +'\n' + \
        str(self.settings['nnls_solver']) + '                                  [ nnls solver ]'

        nn_file= open(path+'ml'+'{:01.2f}'.format(ml)+'/nn.in',"w")
        nn_file.write(text)
        nn_file.close()

    def read_weights(self):
        fname = self.mod_dir_with_ml + '/nn_orb.out'
        col_names = ['orb_idx',
                     'E_idx',
                     'I2_idx',
                     'I3_idx',
                     'totalnotregularizable', # see line 535 of orblib_f.f90
                     'orb_type',
                     'weight',
                     'lcut'] # lines 1321-1322 of triaxnnls_CRcut.f90
        # NOTE: column 'lcut' is not present if different "triaxnnls" file used
        dtype = [int, int, int, int, int, int, np.float64, int]
        weights = np.genfromtxt(fname,
                                skip_header=1,
                                names=col_names,
                                dtype=dtype)
        weights = table.Table(weights)
        self.weights = weights

    def read_nnls_orbmat_noerror(self):
        fname = self.mod_dir_with_ml + '/nn_orbmat_noerror.out'
        orbmat_shape = np.loadtxt(fname, max_rows=1, dtype=int)
        orbmat_noerror = np.loadtxt(fname, skiprows=1)
        orbmat_noerror = np.reshape(orbmat_noerror, orbmat_shape)
        return orbmat_noerror

    def read_nnls_orbmat_rhs_and_solution(self):
        fname = self.mod_dir_with_ml + '/nn_orbmat.out'
        orbmat_shape = np.loadtxt(fname, max_rows=1, dtype=int)
        orbmat_size = np.product(orbmat_shape)
        tmp = np.loadtxt(fname, skiprows=1)
        orbmat = tmp[0:orbmat_size]
        orbmat = np.reshape(orbmat, orbmat_shape)
        rhs = tmp[orbmat_size:orbmat_size+orbmat_shape[1]]
        solution = tmp[orbmat_size+orbmat_shape[1]:]
        return orbmat, rhs, solution

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



class PrashsCoolNewWeightSolver(WeightSolver):

    def __init__(self,
                 system=None,
                 settings=None,
                 directory_noml=None,
                 CRcut=False):
        self.system = system
        self.settings = settings
        self.direc = directory_noml
        self.CRcut = CR_cut
        self.get_observed_constraints()
        self.ennumerate_constraints()

    def get_observed_constraints(self):
        mge = self.system.cmp_list[2].mge
        # intrinsic mass
        intrinsic_masses = mge.get_intrinsic_masses_from_file(self.direc)
        self.intrinsic_masses = intrinsic_masses
        self.intrinsic_mass_error = self.settings['lum_intr_rel_err']
        # projected
        projected_masses = mge.get_projected_masses_from_file(self.direc)
        self.projected_masses = projected_masses
        self.projected_mass_error = self.settings['sb_proj_rel_err']
        # total mass constraint
        self.total_mass = np.sum(intrinsic_masses)
        self.total_mass_error = np.min([self.intrinsic_mass_error/10.,
                                        np.abs(1. - self.total_mass)])
        # get observed gh values and errors
        kinematics = self.system.cmp_list[2].kinematic_data[0]
        tmp = kinematics.get_observed_values_and_uncertainties()
        obs_gh, obs_gh_err = tmp
        # ... and scale both by projected masses
        obs_gh = (obs_gh.T * self.projected_masses).T
        obs_gh_err = (obs_gh_err.T * self.projected_masses).T
        self.obs_gh = obs_gh
        self.obs_gh_err = obs_gh_err

    def ennumerate_constraints(self):
        # enumerate constriants
        n_intrinsic = np.product(self.intrinsic_masses.shape)
        n_apertures = len(self.projected_masses)
        n_gh = self.settings['number_GH']
        # constraints = tot mass, intrinsic masses, projected masses, kinematics
        n_total_constraints = 1 + n_intrinsic + n_apertures + n_gh * n_apertures
        self.n_intrinsic = n_intrinsic
        self.n_apertures = n_apertures
        self.n_gh = n_gh
        self.n_total_constraints = n_total_constraints

    def construct_nnls_matrix_and_rhs(self, orblib):
        # vectors of observed constraints, errors, and orbit properties
        con = np.zeros(self.n_total_constraints)
        econ = np.zeros(self.n_total_constraints)
        orbmat = np.zeros((self.n_total_constraints, orblib.n_orbs))
        # total mass
        con[0] = self.total_mass
        econ[0] = self.total_mass_error
        if econ[0]<=0.0:
            econ[0] = con[0]*0.01
        orbmat[0,:] = 1.
        # intrinsic mass
        idx = slice(1,1+self.n_intrinsic)
        con[idx] = np.ravel(self.intrinsic_masses)
        error = self.intrinsic_masses * self.intrinsic_mass_error
        error = np.abs(np.ravel(error))
        error[np.where(error<=0.)] = 1.0e-16
        econ[idx] = np.abs(np.ravel(error))
        orb_int_masses = orblib.intrinsic_masses
        orb_int_masses = np.reshape(orb_int_masses, (orblib.n_orbs, -1))
        orbmat[idx,:] = orb_int_masses.T
        # projected mass
        idx = slice(1+self.n_intrinsic, 1+self.n_intrinsic+self.n_apertures)
        con[idx] = self.projected_masses
        econ[idx] = np.abs(self.projected_masses * self.projected_mass_error)
        orbmat[idx,:] = orblib.projected_masses.T
        # kinematics
        idx = slice(1+self.n_intrinsic+self.n_apertures, None)
        con[idx] = np.ravel(self.obs_gh.T)
        econ[idx] = np.ravel(self.obs_gh_err.T)
        kinematics = self.system.cmp_list[2].kinematic_data[0]
        # to mimic `triaxnnnls_CRcut.f90`
        # Set the first and last point in the velocity histograms to zero
        orblib.losvd_histograms.y[:,0,:] = 0.
        orblib.losvd_histograms.y[:,-1,:] = 0.
        orb_gh = kinematics.transform_orblib_to_observables(orblib)
        # apply 'CRcut' - cutting orbits where |V - V_obs|> 3sigma_obs
        # see Zhu+2018 MNRAS 2018 473 3000 for details
        if self.CRcut:
            orb_mu_v = orblib.losvd_histograms.get_mean()
            obs_mu_v = kinematics.data['v']
            obs_sig_v = kinematics.data['sigma']
            delta_v = np.abs(orb_mu_v - obs_mu_v)
            condition1 = (np.abs(obs_mu_v)/obs_sig_v > 1.5)
            condition2 = (delta_v/obs_sig_v > 3.0)
            condition3 = (obs_mu_v*orb_mu_v < 0)
            idx_cut = np.where(condition1 & condition2 & condition3)
            cut = np.zeros_like(orb_mu_v, dtype=bool)
            cut[idx_cut] = True
            naperture_cut = np.sum(cut, 1)
            # orbit 'j' is "bad" in naperture_cut[j] apertures
            # if an orbit is bad in 0 or 1 apertures, then we ignore this
            cut[naperture_cut<1,:] = False
            # to cut an orbit, replace it's h1 by 3.0/dvhist(i)
            idx_cut = np.where(cut)
            v_range = float(orblib.settings['hist_vel'])
            n_bins = orblib.losvd_histograms.x.size
            dvhist = v_range/n_bins
            orb_gh[idx_cut[0], idx_cut[1], 0] = 3./dvhist
        orb_gh = np.swapaxes(orb_gh, 1, 2)
        orb_gh = np.reshape(orb_gh, (orblib.n_orbs,-1))
        orbmat[361+152:,:] = orb_gh.T
        # divide constraint vector and matrix by errors
        rhs = con/econ
        orbmat = (orbmat.T/econ).T
        return orbmat, rhs

    def solve(self, orblib):
        A, b = self.construct_nnls_matrix_and_rhs(orblib)
        solution = optimize.nnls(A, b)
        return solution









# end
