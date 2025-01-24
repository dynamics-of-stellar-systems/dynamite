import os
import copy
import numpy as np
from astropy import table
from astropy.io import ascii
import subprocess
import logging
from scipy import optimize

try:
    import cvxopt
except ModuleNotFoundError:
    pass

from dynamite import analysis
from dynamite import physical_system as physys
from dynamite import kinematics as dyn_kin


class WeightSolver(object):
    """Generic WeightSolver class

    Specific implementations are defined as sub-classes. Each one should
    have a main method `solve`

    Parameters
    ----------
    config : a ``dyn.config_reader.Configuration`` object
    directory_with_ml : string
        model directory with the ml extension
    CRcut : Bool, default False
        whether to use the `CRcut` solution for the counter-rotating orbit
        problem. See Zhu et al. 2018 for more. If `CRcut` is given in the
        configuration file's weight solver settings (which is normally the
        case), this parameter is ignored.

    """

    def __init__(self, config, directory_with_ml, CRcut=False):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.config = config
        self.system = config.system
        self.settings = config.settings.weight_solver_settings
        self.direc_with_ml = directory_with_ml
        self.direc_no_ml \
            = directory_with_ml[:directory_with_ml[:-1].rindex('/')+1]
        if 'CRcut' in self.settings.keys():
            CRcut = self.settings['CRcut']
        self.CRcut = CRcut
        self.weight_file = f'{self.direc_with_ml}orbit_weights.ecsv'

    def solve(self, orblib, ignore_existing_weights=False):
        """Template solve method

        Specific implementations should override this.

        Parameters
        ----------
        orblib : dyn.OrbitLibrary object
        ignore_existing_weights : bool
            If True, do not check for already existing weights and solve again.
            Default is False.

        Returns
        -------
        weights : array
            orbit weights
        chi2_all : float
            a total chi2 value
        chi2_kin : float
            a chi2 value purely for kinematics
        chi2_kinmap : float
            directly calculates the chi2 from the kinematic maps
        """
        self.logger.info(f"Using WeightSolver: {__class__.__name__}")
        # ...
        # calculate orbit weights, and model chi2 values here
        # ...
        weights = 0.
        chi2_tot = 0.
        chi2_kin = 0.
        chi2_kinmap = 0.
        # ...
        return weights, chi2_tot, chi2_kin, chi2_kinmap

    def chi2_kinmap(self, weights):
        """
        Returns the chi2 directly calculated from the gh kinematic maps.

        For each kinematic set, the following applies: If number_GH in the
        weight_solver_settings is smaller than the number of GH coefficients
        in the data file, only number_GH coefficients will be considered.
        If number_GH is greater than the number of GH coefficients in the
        data file, only the coefficients in the data file will be considered.

        Does only work with Gauss Hermite kinematics.

        Parameters
        ----------
        weights : ``numpy.array`` like
            The model's orbital weights.

        Returns
        -------
        chi2_kinmap : float
            chi2 directly calculated from the kinematic maps: sum of
            squared residuals of V, sigma, and GH coefficients from h_3 to h_N

        """
        stars = \
          self.system.get_component_from_class(physys.TriaxialVisibleComponent)
        if any(k.type != 'GaussHermite' for k in stars.kinematic_data):
            self.logger.info("'GaussHermite' kinematics required for "
                             "kinmapchi2. Value set to nan.")
            return float('nan')  # #######################################
        number_gh = self.settings['number_GH']
        mod=self.config.all_models.get_model_from_directory(self.direc_with_ml)
        chi2_kinmap = 0.
        for kin_set, kin_data in enumerate(stars.kinematic_data):
            n_gh = min(number_gh, kin_data.max_gh_order)
            coefs = ['v', 'sigma'] + [f'h{i}' for i in range(3, n_gh + 1)]
            # get the model's projected masses=flux (unused) and kinematic data
            a=analysis.Analysis(config=self.config, model=mod, kin_set=kin_set)
            model_gh_coef = a.get_gh_model_kinematic_maps(v_sigma_option='fit',
                                                          weights=weights)
            # get the observed projected masses (unused) and kinematic data
            kinematics_data = \
                kin_data.get_data(self.settings, apply_systematic_error=False)
            # calculate chi2_kinmap
            for coef in coefs:
                obs_val = np.array(kinematics_data[coef])
                mod_val = np.array(model_gh_coef[coef])
                err_val = np.array(kinematics_data['d' + coef])
                chi2_kinmap += sum(np.square((obs_val - mod_val) / err_val))
        return chi2_kinmap

    def weight_file_exists(self):
        """Check whether the file(s) holding the current model's weights exist.

        May be re-implemented by sub-classes.

        Returns
        -------
        bool
            True if weight solving data exists, False otherwise.

        """
        return os.path.isfile(self.weight_file)


class LegacyWeightSolver(WeightSolver):
    """Use `legacy` AKA Fortran weight solving.

    Uses the legcay_fortran program ``triaxnnls_CRcut.f90`` or
    ```triaxnnls_noCRcut.f90``. Uses Lawson and Hanson non-negative
    least-squares algorithm.

    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.legacy_directory=self.config.settings.legacy_settings['directory']
        self.sformat = self.system.parameters[0].sformat # this is ml's format
        ml_idx = self.direc_with_ml.rindex('/ml')
        ml_str = self.direc_with_ml[ml_idx+3:]
        self.ml = float(ml_str[:-1]) if ml_str[-1] == '/' else float(ml_str)
        self.fname_nn_kinem = self.direc_with_ml + 'nn_kinem.out'
        self.fname_nn_nnls = self.direc_with_ml + 'nn_nnls.out'
        # check the format of the orbit library files
        # check == True means there are 2 orblib_* and 2 orblibbox_* files
        # check == False means there are only two orblib files,
        # orblib.dat.bz2 and orbibbox.dat.bz2 (legacy behavior)
        pth = self.direc_no_ml + 'datfil/'
        check = os.path.isfile(f'{pth}orblib_qgrid.dat.bz2') \
                and os.path.isfile(f'{pth}orblib_losvd_hist.dat.bz2') \
                and os.path.isfile(f'{pth}orblibbox_qgrid.dat.bz2') \
                and os.path.isfile(f'{pth}orblibbox_losvd_hist.dat.bz2')
        self.legacy_files = False if check else True
        # prepare fortran input file for nnls
        self.copy_kinematic_data()
        self.create_fortran_input_nnls()

    def copy_kinematic_data(self):
        """Copy kin data to infil/ direc
        """
        stars = self.system.get_unique_triaxial_visible_component()
        kinematics = stars.kinematic_data
        # convert kinematics to old format to input to fortran
        for i in np.arange(len(kinematics)):
            if len(kinematics) == 1:
                old_filename = self.direc_no_ml + 'infil/kin_data.dat'
            else:
                old_filename = self.direc_no_ml+'infil/kin_data_'+str(i)+'.dat'
            kinematics[i].convert_to_old_format(old_filename, self.settings)
        # combine all kinematics into one file
        if len(kinematics) > 1:
            if not all(isinstance(kin, dyn_kin.GaussHermite)
                       for kin in kinematics):
                text = 'Multiple kinematics: all must be GaussHermite'
                self.logger.error(text)
                raise ValueError(text)
            # make a dummy 'kins_combined' object ...
            kins_combined = copy.deepcopy(kinematics[0])
            # ...replace data attribute with stacked table of all kinematics
            kins_combined.data = table.vstack([k.get_data(self.settings)
                                               for k in kinematics])
            kins_combined.n_apertures = len(kins_combined.data)
            kins_combined.max_gh_order = self.settings['number_GH']
            old_filename = self.direc_no_ml + 'infil/kin_data_combined.dat'
            kins_combined.convert_to_old_format(old_filename, self.settings)

    def create_fortran_input_nnls(self):
        """create fortran input file nn.in

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        # When varying ml the LOSVD is scaled - no new orbits are calculated.
        # Therefore we need to know the ml that was used for the orbit library.
        # The scaling factor is sqrt(model_ml/original_orblib_ml).
        mod=self.config.all_models.get_model_from_directory(self.direc_with_ml)
        ml_scaling_factor = \
            self.config.all_models.get_model_velocity_scaling_factor(model=mod)
        #-------------------
        #write nn.in
        #-------------------
        n_kin = len(self.system.get_unique_triaxial_visible_component().kinematic_data)

        if n_kin==1:
            kin_data_file='kin_data.dat'

        else:
            kin_data_file='kin_data_combined.dat'

        text='infil/parameters_pot.in' +'\n' + \
        str(self.settings['regularisation'])   + '                                  [ regularization strength, 0 = no regularization ]' +'\n'  + \
        f'ml{self.ml:{self.sformat}}/nn\n' + \
        'datfil/mass_qgrid.dat' +'\n' + \
        'datfil/mass_aper.dat' +'\n' + \
        str(self.settings['number_GH']) + '	                           [ # of GH moments to constrain the model]' +'\n' + \
        'infil/'+kin_data_file+'\n' + \
        str(self.settings['GH_sys_err']) + '    [ systemic error of v, sigma, h3, h4... ]' + '\n' + \
        str(self.settings['lum_intr_rel_err']) + '                               [ relative error for intrinsic luminosity ]' +'\n' + \
        str(self.settings['sb_proj_rel_err']) + '                               [ relative error for projected SB ]' + '\n' + \
        str(ml_scaling_factor)  + '                                [ scale factor related to M/L, sqrt( (M/L)_k / (M/L)_ref ) ]' + '\n'
        if self.legacy_files:
            text += 2 * f'datfil/orblib_{self.ml}.dat\n' + \
                    2 * f'datfil/orblibbox_{self.ml}.dat\n'  # yes, really...
        else:
            for f in '_qgrid', '_losvd_hist', 'box_qgrid', 'box_losvd_hist':
                text += f'datfil/orblib{f}_{self.ml}.dat\n'
        text += str(self.settings['nnls_solver']) + '                                  [ nnls solver ]'

        nn_file = open(self.direc_no_ml + f'ml{self.ml:{self.sformat}}/nn.in',
                       "w")
        nn_file.write(text)
        nn_file.close()

    def solve(self, orblib=None, ignore_existing_weights=False):
        """Main method to solve NNLS problem.

        Parameters
        ----------
        orblib : dyn.OrbitLibrary
            This parameter is not used in this Legacy implementation (as all
            orbit library information is read from files). It is included here
            for consistency with later WeightSolver implementations
        ignore_existing_weights : bool
            If True, do not check for already existing weights and solve again.
            Default is False.

        Returns
        -------
        tuple
            (weights, chi2_all, chi2_kin, chi2_kinmap) where:
                -   weights : array, of orbit weights
                -   chi2_all : float, sum of squared residuals for intrinsic
                    masses, projected_masses and GH coefficients from h_1 to h_n
                -   chi2_kin : float sum of squared residuals for GH
                    coefficients h_1 to h_n
                -   chi2_kinmap : directly calculates the chi2 from the
                    kinematic maps

        """
        self.logger.info(f"Using WeightSolver: {__class__.__name__}")
        if (not ignore_existing_weights) and self.weight_file_exists():
            self.logger.info("Reading NNLS solution from existing output "
                             f"{self.weight_file}.")
            results = ascii.read(self.weight_file)
            weights = results['weights']
            chi2_tot = results.meta['chi2_tot']
            chi2_kin = results.meta['chi2_kin']
            chi2_kinmap = results.meta['chi2_kinmap']
        else:
            fname_nn_orbmat = self.direc_with_ml + 'nn_orbmat.out'
            # If legacy result files do not exist, run weight solving.
            check = (os.path.isfile(self.fname_nn_kinem) and
                     os.path.isfile(self.fname_nn_nnls) and
                     os.path.isfile(fname_nn_orbmat))
            if ignore_existing_weights or not check:
                for f in [self.fname_nn_kinem, self.fname_nn_nnls,
                          fname_nn_orbmat]:
                    if os.path.isfile(f):
                        os.remove(f)
                # set the current directory to the directory in which
                # the models are computed
                cur_dir = os.getcwd()
                os.chdir(self.direc_no_ml)
                cmdstr = self.write_executable_for_weight_solver()
                with open(cmdstr) as f:
                    for line in f:
                        i = line.find('>>')
                        if i >= 0:
                            j = line.find('.log')
                            logfile = line[i+3:j+4]
                            break
                self.logger.info("Fitting orbit library to the kinematic " +
                                 f"data: {logfile[:logfile.rindex('/')]}")
                p = subprocess.run('bash '+cmdstr,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT,
                                   shell=True)
                # clean up decompressed files
                if self.legacy_files:
                    for f_name in [f'datfil/orblib_{self.ml}.dat',
                                f'datfil/orblibbox_{self.ml}.dat']:
                        if os.path.isfile(f_name):
                            os.remove(f_name)
                else:
                    for f in '_qgrid', '_losvd_hist', 'box_qgrid', 'box_losvd_hist':
                        f_name = f'datfil/orblib{f}_{self.ml}.dat'
                        if os.path.isfile(f_name):
                            os.remove(f_name)
                log_file = f'Logfile: {self.direc_no_ml+logfile}.'
                if not p.stdout.decode("UTF-8"):
                    self.logger.info(f'...done, NNLS problem solved - {cmdstr}'
                                     f' exit code {p.returncode}. {log_file}')
                else:
                    text = f'...failed! {cmdstr} exit code {p.returncode}. ' \
                           f'Message: {p.stdout.decode("UTF-8")}'
                    if p.returncode == 127: # command not found
                        text += 'Check DYNAMITE legacy_fortran executables.'
                        self.logger.error(text)
                        raise FileNotFoundError(text)
                    text += f'{log_file} Be wary: DYNAMITE may crash...'
                    self.logger.warning(text)
                    raise RuntimeError(text)
                # set the current directory to the dynamite directory
                os.chdir(cur_dir)
                # delete existing .yaml files and copy current config file
                # into model directory
                self.config.copy_config_file(self.direc_with_ml)
            else:  # If legacy output files exist, just create the weight file
                self.logger.info("Reading NNLS solution from existing legacy "
                                 "output and converting to weights file.")
            # Now the legacy result files exist -> read, calculate
            # kinmapchi2, and save to the weight file.
            weights, chi2_tot, chi2_kin = \
                self.get_weights_and_chi2_from_orbmat_file()
            chi2_kinmap = self.chi2_kinmap(weights)
            results = table.Table()
            results['weights'] = weights
            results.meta = {'chi2_tot': chi2_tot,
                            'chi2_kin': chi2_kin,
                            'chi2_kinmap': chi2_kinmap}
            results.write(self.weight_file,
                          format='ascii.ecsv',
                          overwrite=True)
            # clean up
            if os.path.isfile(fname_nn_orbmat):
                os.remove(fname_nn_orbmat)
        return weights, chi2_tot, chi2_kin, chi2_kinmap

    def write_executable_for_weight_solver(self):
        """write executable bash script file

        Parameters
        ----------
        None

        Returns
        -------
        string
            the name of the bash script file to execute

        """
        nn = f'ml{self.ml:{self.sformat}}/nn'
        cmdstr = f'cmd_nnls_{self.ml}'
        txt_file = open(cmdstr, "w")
        txt_file.write('#!/bin/bash' + '\n')
        txt_file.write('# if the gzipped orbit library exist unzip it' + '\n')
        if self.legacy_files:
            txt_file.write(f'test -e datfil/orblib_{self.ml}.dat || bunzip2 -c  datfil/orblib.dat.bz2 > datfil/orblib_{self.ml}.dat' + '\n')
            txt_file.write(f'test -e datfil/orblibbox_{self.ml}.dat || bunzip2 -c  datfil/orblibbox.dat.bz2 > datfil/orblibbox_{self.ml}.dat' + '\n')
        else:
            for f in '_qgrid', '_losvd_hist', 'box_qgrid', 'box_losvd_hist':
                file_name = f'datfil/orblib{f}_{self.ml}.dat'
                txt_file.write(f'test -e {file_name} || '
                               f'bunzip2 -c  datfil/orblib{f}.dat.bz2 > {file_name}\n')
        if self.system.is_bar_disk_system():
            txt_file.write(f'test -e {self.legacy_directory}/triaxnnls_bar' +
                           f' || {{ echo "File {self.legacy_directory}/triaxnnls_bar not found." && exit 127; }}\n')
            txt_file.write('test -e ' + str(nn) + '_kinem.out || ' +
                           self.legacy_directory +
                           f'/triaxnnls_bar < {nn}.in >> {nn}ls.log '
                           '|| exit 1\n')
        elif self.CRcut is True:
            txt_file.write(f'test -e {self.legacy_directory}/triaxnnls_CRcut' +
                           f' || {{ echo "File {self.legacy_directory}/triaxnnls_CRcut not found." && exit 127; }}\n')
            txt_file.write('test -e ' + str(nn) + '_kinem.out || ' +
                           self.legacy_directory +
                           f'/triaxnnls_CRcut < {nn}.in >> {nn}ls.log '
                           '|| exit 1\n')
        else:
            txt_file.write(f'test -e {self.legacy_directory}/triaxnnls_noCRcut' +
                           f' || {{ echo "File {self.legacy_directory}/triaxnnls_noCRcut not found." && exit 127; }}\n')
            txt_file.write('test -e ' + str(nn) + '_kinem.out || ' +
                           self.legacy_directory +
                           f'/triaxnnls_noCRcut < {nn}.in >> {nn}ls.log '
                           '|| exit 1\n')
        txt_file.close()
        return cmdstr

    def read_weights(self):
        """Read ``nn_orb.out`` to astropy table

        this contains oribtal weights, orbit type, and other columns

        Returns
        -------
        None
            sets ``self.weights`` which is an astropy table containing
            the orbital weights

        """
        fname = self.direc_with_ml + 'nn_orb.out'
        col_names = ['orb_idx',
                     'E_idx',
                     'I2_idx',
                     'I3_idx',
                     'totalnotregularizable', # see line 535 of orblib_f.f90
                     'orb_type',
                     'weight',
                     'lcut'] # lines 1321-1322 of triaxnnls_CRcut.f90
        # NOTE: column 'lcut' is not present if different "triaxnnls" file used
        dtype = [int, int, int, int, int, int, float, int]
        weights = np.genfromtxt(fname,
                                skip_header=1,
                                names=col_names,
                                dtype=dtype)
        weights = table.Table(weights)
        self.weights = weights

    def read_nnls_orbmat_rhs_and_solution(self):
        """Read ``nn_orbmat.out``

        This contains the matrix and right-hand-side for the NNLS problem, and
        the solution

        Returns
        -------
        tuple
            (orbmat, rhs, solution)

        """
        fname = self.direc_with_ml + 'nn_orbmat.out'
        orbmat_shape = np.loadtxt(fname, max_rows=1, dtype=int)
        orbmat_size = np.prod(orbmat_shape)
        tmp = np.loadtxt(fname, skiprows=1)
        orbmat = tmp[0:orbmat_size]
        orbmat = np.reshape(orbmat, orbmat_shape)
        orbmat = orbmat.T
        rhs = tmp[orbmat_size:orbmat_size+orbmat_shape[1]]
        solution = tmp[orbmat_size+orbmat_shape[1]:]
        return orbmat, rhs, solution

    def get_weights_and_chi2_from_orbmat_file(self):
        """
        Get weights and chi2 from ``nn_orbmat.out``

        **Note**: Chi2 values returned differ from `read_chi2` method.
        See that docstring for more.

        Returns
        -------
        tuple
            (weights, chi2_all, chi2_gh), where:

                -   weights : array of orbit weights
                -   chi2_all : sum of squared residuals for intrinsic masses,
                    projected_masses and GH coefficients h_1 to h_n
                -   chi2_kin : sum of squared residuals for GH coefficients h_1 to h_n

        """
        A, b, weights = self.read_nnls_orbmat_rhs_and_solution()
        chi2_vector = (np.dot(A, weights) - b)**2.
        chi2_tot = np.sum(chi2_vector)

        if self.system.is_bar_disk_system():
            bardisk = self.system.get_unique_bar_component()
            mge = bardisk.mge_lum + bardisk.disk_lum
        else:
            stars = self.system.get_unique_triaxial_visible_component()
            mge = stars.mge_lum

        intrinsic_masses = mge.get_intrinsic_masses_from_file(self.direc_no_ml)
        projected_masses = mge.get_projected_masses_from_file(self.direc_no_ml)
        n_intrinsic = np.prod(intrinsic_masses.shape)
        n_apertures = len(projected_masses)
        chi2_kin = np.sum(chi2_vector[1+n_intrinsic+n_apertures:])
        return weights, chi2_tot, chi2_kin

    # def chi2_kinmap(self):
    #     """
    #     Returns the chi2 directly calculated from the kinematic maps.

    #     Returns
    #     -------
    #     chi2_kinmap : float
    #         chi2 directly calculated from the kinematic maps.

    #     """
    #     _, chi2_kinmap = self.read_chi2()
    #     return chi2_kinmap

    def read_chi2(self):
        """Read chi2 values from `nn_kinem.out`

        Taken from old `schwpy` code, lines 181-212 of schw_domoditer.py

        **Note**:
        This is a legacy method for reading legacy output and it not used by
        default. Instead we use ``self.get_chi2_from_orbmat`` get chi2 values.
        The chi2 value definitions of this method are NOT the same chi2 values
        given by ``self.get_chi2_from_orbmat``. They differ in
        (i) including intrinsic/projected mass constraints, and (ii) using
        h1/h2 vs V/sigma, and (iii) if CRcut==True, whether the 'cut' orbits
        - with artificially large h1 - are included (here they aren't)

        Returns
        -------
        tuple
            (chi2, kinchi2) where:
                -   chi2 = sum of sq. residuals of observed GH coefficients h_1
                    to h_N
                -   kinchi2 = sum of sq. residuals of V, sigma, and GH
                    coefficients from h_3 to h_N

        """
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
        # k = is array of column indices, for [V, sigma, h3, ..., h_ngh]
        #                       observed   modelled        error
        kinchi2 = sum(sum(pow(((ka[:, k] - ka[:, k + 1]) / ka[:, k + 2]), 2.0)))
        return chi2, kinchi2

    def __read_file_element(self, infile, rows, cols):
        """Read fields in a tabular data according to the their row/column.

        Taken from schwpy schw_misc

        Parameters
        ----------
        infile : string
            input file
        rows : array of ints
            row array of locations indexed starts from 1
        cols : array of ints
            column array of locations indexed starts from 1

        Returns
        -------
        array read from given locations in file

        """
        lines = [line.rstrip('\n').split() for line in open(infile)]
        output=[]
        for i in range(0, len(rows)):
            output.append(lines[rows[i] - 1][cols[i] - 1])
        return output


class NNLS(WeightSolver):
    """Python implementations of NNLS weight solving

    Uses either scipy.optimize.nnls or cvxopt as backends. This constructs the
    NNLS matrix and rhs, solves, and saves the result.

    Parameters
    ----------
    nnls_solver : string
        either ``scipy`` or ``cvxopt``

    """
    def __init__(self,
                 nnls_solver=None,
                 **kwargs):
        super().__init__(**kwargs)
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        if nnls_solver is None:
            nnls_solver = self.settings['nnls_solver']
        assert nnls_solver in ['scipy', 'cvxopt'], 'Unknown nnls_solver'
        self.nnls_solver = nnls_solver
        self.get_observed_mass_constraints()

    def get_observed_mass_constraints(self):
        """Get aperture+intrinsic mass constraints from MGE

        Returns the projected masses of the MGE for the kinematic data
        apertures.

        Returns
        -------
        None
            sets attributes:

                - ``self.intrinsic_masses``
                - ``self.intrinsic_mass_error``
                - ``self.projected_masses``
                - ``self.projected_mass_error``
                -   constraint counts ``self.n_intrinsic``, ``self.n_apertures``
                    and ``self.n_mass_constraints``

        """
        if self.system.is_bar_disk_system():
            bardisk = self.system.get_unique_bar_component()
            mge = bardisk.mge_lum + bardisk.disk_lum
        else:
            stars = self.system.get_unique_triaxial_visible_component()
            mge = stars.mge_lum

        # intrinsic mass
        intrinsic_masses = mge.get_intrinsic_masses_from_file(self.direc_no_ml)
        self.intrinsic_masses = intrinsic_masses
        self.intrinsic_mass_error = self.settings['lum_intr_rel_err']
        # projected
        projected_masses = mge.get_projected_masses_from_file(self.direc_no_ml)
        self.projected_masses = projected_masses
        self.projected_mass_error = self.settings['sb_proj_rel_err']
        # total mass constraint
        self.total_mass = np.sum(intrinsic_masses)
        self.total_mass_error = np.min([self.intrinsic_mass_error/10.,
                                        np.abs(1. - self.total_mass)])
        # enumerate the mass constriants
        n_intrinsic = np.prod(self.intrinsic_masses.shape)
        n_apertures = len(self.projected_masses)
        self.n_intrinsic = n_intrinsic
        self.n_apertures = n_apertures
        # mass constraints = total mass (1) + intrinsic mass + aperture mass
        self.n_mass_constraints = 1 + n_intrinsic + n_apertures

    def construct_nnls_matrix_and_rhs(self, orblib):
        """construct nnls matrix_and rhs

        Parameters
        ----------
        orblib : ``dyn.orblib.OrbitLibrary``
            an orbit library

        Returns
        -------
        tuple
            (orbmat, rhs)

        """
        # construct vector of observed constraints (con), errors (econ) and
        # matrix or orbit propertites (orbmat)
        con = np.zeros(self.n_mass_constraints)
        econ = np.zeros(self.n_mass_constraints)
        orbmat = np.zeros((self.n_mass_constraints, orblib.n_orbs))
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
        orbmat[idx,:] = np.hstack(orblib.projected_masses).T
        # add kinematics to con, econ, orbmat
        stars = self.system.get_unique_triaxial_visible_component()
        kins_and_orb_losvds = zip(stars.kinematic_data, orblib.losvd_histograms)
        idx_ap_start = 0
        for (kins, orb_losvd) in kins_and_orb_losvds:
            # pick out the projected masses for this kinematic set
            n_ap = kins.n_spatial_bins  # OK for both GaussHermite & BayesLOSVD
            idx_ap_end = idx_ap_start + n_ap
            prj_mass_i = self.projected_masses[idx_ap_start:idx_ap_end]
            idx_ap_start += n_ap
            # scale observed kinematics and errors by projected masses
            tmp = kins.get_observed_values_and_uncertainties(self.settings)
            obs_kins, obs_kins_err = tmp
            obs_kins = (obs_kins.T * prj_mass_i).T
            obs_kins_err = (obs_kins_err.T * prj_mass_i).T
            # set the first and last point in the velocity histograms to zero
            # to mimic what is done in `triaxnnnls_CRcut.f90`
            orb_losvd.y[:,0,:] = 0.
            orb_losvd.y[:,-1,:] = 0.
            # transform orblib to same parameterisation as observed kinematics
            orb_kins = kins.transform_orblib_to_observables(orb_losvd,
                                                            self.settings)
            if self.CRcut:
                # note: this only has an effect if type(kins) is GaussHermite
                orb_kins = self.apply_CR_cut(kins, orb_losvd, orb_kins)
            # append constraints/errors/orbits to con/econ/orbmat
            obs_kins = np.ravel(obs_kins)
            con = np.concatenate((con, obs_kins))
            obs_kins_err = np.ravel(obs_kins_err)
            econ = np.concatenate((econ, obs_kins_err))
            orb_kins = np.reshape(orb_kins, (orblib.n_orbs, -1))
            orbmat = np.vstack((orbmat, orb_kins.T))
        # divide constraint vector and matrix by errors
        rhs = con/econ
        orbmat = (orbmat.T/econ).T
        return orbmat, rhs

    def apply_CR_cut(self, kins, orb_losvd, orb_gh):
        r"""apply `CRcut`

        to solve the `counter rotating orbit problem`. This cuts orbits which
        have :math:`|V - V_\mathrm{obs}|> 3\sigma_\mathrm{obs}`. See
        Zhu+2018 MNRAS 2018 473 3000 for details

        Parameters
        ----------
        kins : a ``dyn.kinematics.Kinematic`` object
        orb_losvd : ``dyn.kinematics.Histogram``
            historgram of orblib losvds
        orb_gh : array
            array of input gh expansion coefficients, before the CRcut

        Returns
        -------
        array
            array of input gh expansion coefficients, after the CRcut

        """
        if type(kins) is not dyn_kin.GaussHermite:
            return orb_gh
        orb_mu_v = orb_losvd.get_mean()
        kins_data = kins.get_data(self.settings, apply_systematic_error=False)
        obs_mu_v = kins_data['v']
        obs_sig_v = kins_data['sigma']
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
        dvhist = kins.hist_width/kins.hist_bins
        dvhist = np.max(dvhist)
        orb_gh[idx_cut[0], idx_cut[1], 0] = 3./dvhist
        return orb_gh

    def solve(self, orblib, ignore_existing_weights=False):
        """Solve for orbit weights

        **Note:** the returned chi2 values are not the same as
        ``LegacyWeightSolver.read_chi2`` - see the docstring for more info

        Parameters
        ----------
        orblib : dyn.OrbitLibrary
        ignore_existing_weights : bool
            If True, do not check for already existing weights and solve again.
            Default is False.

        Returns
        -------
        tuple
            (weights, chi2_all, chi2_kin, chi2_kinmap) where:
                -   weights : array, of orbit weights
                -   chi2_all : float, sum of squared residuals for intrinsic
                    masses, projected_masses and GH coefficients from h_1 to h_n
                -   chi2_kin : float sum of squared residuals for GH
                    coefficients h_1 to h_n
                -   chi2_kinmap : directly calculates the chi2 from the
                    kinematic maps

        """
        self.logger.info(f"Using WeightSolver: {__class__.__name__}/"
                         f"{self.nnls_solver}")
        orblib.read_losvd_histograms()  # sets orblib.losvd_histograms,
                                        # orblib.intrinsic_masses, and
                                        # orblib.projected_masses
        if (not ignore_existing_weights) and self.weight_file_exists():
            results = ascii.read(self.weight_file, format='ecsv')
            self.logger.info("NNLS solution read from existing output "
                             f"{self.weight_file}.")
            weights = results['weights']
            chi2_tot = results.meta['chi2_tot']
            chi2_kin = results.meta['chi2_kin']
            chi2_kinmap = results.meta['chi2_kinmap']
        else:
            A, b = self.construct_nnls_matrix_and_rhs(orblib)
            if self.nnls_solver=='scipy':
                try:
                    solution = optimize.nnls(A, b)
                    weights = solution[0]
                except Exception as e:
                    txt = f'Orblib {orblib.mod_dir}, ml={orblib.parset["ml"]}'\
                        f': SciPy solver error occured: {e} All weights ' \
                        'and chi2 set to nan. Consider trying cvxopt.'
                    self.logger.warning(txt)
                    weights = np.full(A.shape[1], np.nan)
            elif self.nnls_solver=='cvxopt':
                try:
                    P = np.dot(A.T, A)
                    q = -1.*np.dot(A.T, b)
                    solver = CvxoptNonNegSolver(P, q)
                    weights = solver.beta
                except Exception as e:
                    txt = f'Orblib {orblib.mod_dir}, ml={orblib.parset["ml"]}'\
                        f': CVXOPT solver error occured: {e} All weights ' \
                        'and chi2 set to nan. Consider trying scipy.'
                    self.logger.warning(txt)
                    weights = np.full(A.shape[1], np.nan)
            else:
                text = 'Unknown nnls_solver'
                self.logger.error(text)
                raise ValueError(text)
            if not np.isnan(weights[0]):
                # calculate chi2s
                chi2_vector = (np.dot(A, weights) - b)**2.
                chi2_tot = np.sum(chi2_vector)
                chi2_kin = np.sum(chi2_vector[self.n_mass_constraints:])
                chi2_kinmap = self.chi2_kinmap(weights)
                # save the output
                results = table.Table()
                results['weights'] = weights
                # add chi2 to meta data
                results.meta = {'chi2_tot': chi2_tot,
                                'chi2_kin': chi2_kin,
                                'chi2_kinmap': chi2_kinmap}
                results.write(self.weight_file,
                              format='ascii.ecsv',
                              overwrite=True)
                self.logger.info("NNLS problem solved and chi2 calculated for "
                                 f"model {self.direc_with_ml}.")
            else:
                chi2_tot = chi2_kin = chi2_kinmap = np.nan
            # delete existing .yaml files and copy current config file
            # into model directory
            self.config.copy_config_file(self.direc_with_ml)
        return weights, chi2_tot, chi2_kin, chi2_kinmap


class CvxoptNonNegSolver():
    """Solver for NNLS problem using CVXOPT

    Solves the QP problem:
        argmin (1/2 beta^T P beta + q beta T)
        subject to (component-wise) beta > 0

    Parameters
    ----------
    P : array (p, p)
        quadratic part of objective function
    q : array (p,)
        linear part of objective function

    Attributes
    ----------
    success : bool
        whether solver was successful
    beta : array (p,)
        solution

    """

    def __init__(self, P=None, q=None):
        p = P.shape[0]
        P = cvxopt.matrix(P)
        q = cvxopt.matrix(q)
        G = cvxopt.matrix(-np.identity(p))
        h = cvxopt.matrix(np.zeros(p))
        sol = cvxopt.solvers.qp(P, q, G, h)
        self.success = sol['status']=='optimal'
        self.beta = np.squeeze(np.array(sol['x']))

# end
