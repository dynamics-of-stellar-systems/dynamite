import os
import subprocess
import shutil
import logging
import numpy as np
from scipy.io import FortranFile
from astropy import table
import astropy.units as u
import matplotlib.pyplot as plt
import sparse

from dynamite import physical_system as physys
from dynamite import kinematics as dyn_kin

class OrbitLibrary(object):
    """An abstract class for orbit libraries.

    Parameters
    ----------
    system : a ``dyn.physical_system.System`` object
    settings : a ``dyn.config_reader.settings`` object

    """
    def __init__(self,
                 config=None):
        self.system = config.system
        self.settings = config.settings.orblib_settings
        self.generate_ics()
        # self.integrate_loop(timesteps)

    def generate_ics(self):
        #ics: initial conditions
        pass

    def integrate_orbits(self):
        pass


class LegacyOrbitLibrary(OrbitLibrary):
    """Orbit libraries calculated from `legacy` Fortan programs

    """
    def __init__(self,
                 config=None,
                 mod_dir=None,
                 parset=None):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        if config is None:
            text = f'{__class__.__name__} needs configuration object, ' \
                   'None provided.'
            self.logger.error(text)
            raise ValueError(text)
        self.mod_dir = mod_dir
        self.parset = parset
        self.system = config.system
        self.settings = config.settings.orblib_settings
        self.legacy_directory = config.settings.legacy_settings['directory']
        self.in_dir = config.settings.io_settings['input_directory']
        self.orblibs_in_parallel = \
            config.settings.multiprocessing_settings['orblibs_in_parallel']

    def get_orblib(self):
        """main method to calculate orbit libraries

        Writes and executes bash scripts to (i) calculate orbit initial
        conditions, (ii) calculate orbit libraries, (iii) calculate aperture and
        3D grid masses for the MGE. If orbit libraries for this model already
        exist, then this method does nothing.

        Returns
        -------
        Creates the following output files in ``output/models/*/datfil/``:
            - begin.dat                     (ics for tube orbits)
            - beginbox.dat                  (ics for box orbits)
            - orblib.dat.bz2                (zipped tube orbit library)
            - orblib.dat_orbclass.out       (orbit classification for tube orbs)
            - orblibbox.dat.bz2             (zipped box orbit library)
            - orblibbox.dat_orbclass.out    (orbit classification for box orbs)
            - mass_aper.dat                 (MGE masses in apertures)
            - mass_qgrid.dat                (MGE masses in 3D grid)
            - mass_radmass.dat              (MGE masses in radial bins)
            - +8 log and status files

        """
        # check if orbit library was calculated already
        check1 = os.path.isfile(self.mod_dir+'datfil/orblib.dat.bz2')
        check2 = os.path.isfile(self.mod_dir+'datfil/orblibbox.dat.bz2')
        if not check1 or not check2:
            # prepare the fortran input files for orblib
            self.create_fortran_input_orblib(self.mod_dir+'infil/')
            stars = self.system.get_component_from_class( \
                                            physys.TriaxialVisibleComponent)
            kinematics = stars.kinematic_data
            # create the kinematic input files for each kinematic dataset
            for i in np.arange(len(kinematics)):
                # copy aperture and bins file across
                aperture_file = self.in_dir + kinematics[i].aperturefile
                shutil.copyfile(aperture_file,
                            self.mod_dir+'infil/'+ kinematics[i].aperturefile)
                binfile = self.in_dir + kinematics[i].binfile
                shutil.copyfile(binfile,
                            self.mod_dir+'infil/'+ kinematics[i].binfile)
            # calculate orbit libary
            self.get_orbit_ics()
            if self.orblibs_in_parallel:
                self.get_orbit_library_par()
            else:
                self.get_orbit_library()

    def create_fortran_input_orblib(self, path):
        """write input files for Fortran orbit library programs

        Parameters
        ----------
        path : string
            path of the model's ``infil`` directory

        Returns
        -------
        Cretaes the following files in the ``infil`` directory:
            - parameters_pot.in
            - parameters_lum.in
            - orblib.in
            - orblibbox.in
            - triaxmass.in
            - triaxmassbin.in

        """
        #---------------------------------------------
        #write parameters_pot.in and parameters_lum.in
        #---------------------------------------------
        stars = \
          self.system.get_component_from_class(physys.TriaxialVisibleComponent)
        bh = self.system.get_component_from_class(physys.Plummer)
        # used to derive the viewing angles
        q = self.parset[f'q-{stars.name}']
        p = self.parset[f'p-{stars.name}']
        u = self.parset[f'u-{stars.name}']
        theta, psi, phi = stars.triax_pqu2tpp(p,q,u)
        # get dark halo
        dh = self.system.get_all_dark_non_plummer_components()
        self.logger.debug('Checking number of non-plummer dark components')
        error_msg = 'only one non-plummer dark component should be present'
        assert len(dh)==1, error_msg
        self.logger.debug('...checks ok.')
        dh = dh[0]  # extract the one and only dm component

        if isinstance(dh, physys.NFW_m200_c):
            #fix c via m200_c relation, for legacy Fortran it is still NFW
            dm_specs, dm_par_vals = dh.get_dh_legacy_strings(self.parset,
                                                             self.system)
        else:
            dm_specs, dm_par_vals = dh.get_dh_legacy_strings(self.parset)

        # header
        len_mge_pot = len(stars.mge_pot.data)
        len_mge_lum = len(stars.mge_lum.data)
        settngs = self.settings
        text = f'{self.system.distMPc}\n'
        text += f'{theta:06.9f} {phi:06.9f} {psi:06.9f}\n'
        text += f"{self.parset['ml']}\n"
        text += f"{self.parset[f'm-{bh.name}']}\n"
        text += f"{self.parset[f'a-{bh.name}']}\n"
        text += f"{settngs['nE']} {settngs['logrmin']} {settngs['logrmax']}\n"
        text += f"{settngs['nI2']}\n"
        text += f"{settngs['nI3']}\n"
        text += f"{settngs['dithering']}\n"
        text += f"{dm_specs}\n"
        text += f"{dm_par_vals}"
        # parameters_pot.in
        np.savetxt(path + 'parameters_pot.in',
                   stars.mge_pot.data,
                   header=str(len_mge_pot),
                   footer=text,
                   comments='',
                   fmt=['%10.2f','%10.5f','%10.5f','%10.2f'])
        # parameters_lum.in
        np.savetxt(path + 'parameters_lum.in',
                   stars.mge_lum.data,
                   header=str(len_mge_lum),
                   footer=text,
                   comments='',
                   fmt=['%10.2f','%10.5f','%10.5f','%10.2f'])
        #-------------------
        #write orbstart.in
        #-------------------
        text = f"{self.settings['random_seed']}\n"
        text += 'infil/parameters_pot.in' +'\n' + \
        'datfil/orbstart.dat' +'\n' + \
        'datfil/begin.dat' +'\n' + \
        'datfil/beginbox.dat'
        orbstart_file= open(path+'orbstart.in',"w")
        orbstart_file.write(text)
        orbstart_file.close()
        #---------------------------------
        #write orblib.in and orblibbox.in
        #---------------------------------
        def write_orblib_dot_in(box=False):
            tab = '\t\t\t\t\t\t\t\t'
            if box:
                f = open(path +'orblibbox.in', 'w')
            else:
                f = open(path +'orblib.in', 'w')
            f.write(f"{self.settings['random_seed']}\n")
            f.write('#counterrotation_setupfile_version_1\n')
            f.write('infil/parameters_pot.in\n')
            if box:
                f.write(f'datfil/beginbox.dat\n')
            else:
                f.write(f'datfil/begin.dat\n')
            label = '[number of orbital periods to integrate]'
            line = f"{self.settings['orbital_periods']}\t\t{label}\n"
            f.write(line)
            label = '[points to sample for each orbit in the merid. plane]'
            line = f"{self.settings['sampling']}{tab}{label}\n"
            f.write(line)
            label = '[starting orbit]'
            line = f"{self.settings['starting_orbit']}{tab}{label}\n"
            f.write(line)
            label = '[orbits  to intergrate; -1 --> all orbits]'
            line = f"{self.settings['number_orbits']}{tab}{label}\n"
            f.write(line)
            label = '[accuracy]'
            line = f"{self.settings['accuracy']}{tab}{label}\n"
            f.write(line)
            n_psf = len(stars.kinematic_data)
            label = '[number of psfs of the kinematic data]'
            line = f"{n_psf}{tab}{label}\n"
            f.write(line)
            for i in range(n_psf):
                psf_i = stars.kinematic_data[i].PSF
                n_gauss_psf_i = len(psf_i['sigma'])
                label = f'[# of gaussians in psf {i+1}]'
                line = f"{n_gauss_psf_i}{tab}{label}\n"
                f.write(line)
            for i in range(n_psf):
                psf_i = stars.kinematic_data[i].PSF
                for j in range(n_gauss_psf_i):
                    weight_ij, sigma_ij = psf_i['weight'][j], psf_i['sigma'][j]
                    label = f'[weight, sigma of comp {j+1} of psf {i+1}]'
                    line = f"{weight_ij} {sigma_ij}{tab}{label}\n"
                    f.write(line)
            # every kinematic-set has a psf and aperture, so n_psf = n_aperture
            n_aperture = n_psf
            label = '[# of apertures]'
            line = f'{n_aperture}{tab}{label}\n'
            f.write(line)
            for i in np.arange(n_aperture):
                kin_i = stars.kinematic_data[i]
                # aperturefile cant have a label since fortran reads whole line
                f.write(f'"infil/{kin_i.aperturefile}"\n')
                label = f'[use psf {i+1} for {kin_i.name}]'
                line = f'{i+1}{tab}{label}\n'
                f.write(line)
            for i in np.arange(n_aperture):
                kin_i = stars.kinematic_data[i]
                label = f'[vhist width, center and nbins for {kin_i.name}]'
                w, c, b = kin_i.hist_width, kin_i.hist_center, kin_i.hist_bins
                line = f'{w} {c} {b}{tab}{label}\n'
                f.write(line)
            for i in np.arange(n_aperture):
                label = f'[use binning for aperture {1+i}? 0/1 = yes/no]'
                line = f'1{tab}{label}\n'
                f.write(line)
            for i in np.arange(n_aperture):
                kin_i = stars.kinematic_data[i]
                label = f'[binfile for aperture {1+i}]'
                line = f'"infil/{kin_i.binfile}"{tab}{label}\n'
                f.write(line)
            if box:
                f.write(f'datfil/orblibbox.dat\n')
            else:
                f.write(f'datfil/orblib.dat\n')
            f.close()
        write_orblib_dot_in(box=False)
        write_orblib_dot_in(box=True)
        #-------------------
        #write triaxmass.in
        #-------------------
        text='infil/parameters_lum.in' +'\n' + \
        'datfil/orblib.dat' +'\n' + \
        'datfil/mass_radmass.dat' +'\n' + \
        'datfil/mass_qgrid.dat'
        triaxmass_file= open(path+'triaxmass.in',"w")
        triaxmass_file.write(text)
        triaxmass_file.close()
        #-----------------------
        #write triaxmassbin.in
        #-----------------------
        tab = '\t\t\t\t\t\t\t\t'
        n_psf = len(stars.kinematic_data)
        f = open(path + 'triaxmassbin.in', 'w')
        f.write('infil/parameters_lum.in\n')
        f.write(f'{n_psf}{tab}[# of apertures]\n')
        for i in range(n_psf):
            kin_i = stars.kinematic_data[i]
            f.write(f'"infil/{kin_i.aperturefile}"\n')
            psf_i = kin_i.PSF
            n_gauss_psf_i = len(psf_i['sigma'])
            label = f'[# of gaussians in psf {i+1}]'
            line = f"{n_gauss_psf_i}{tab}{label}\n"
            f.write(line)
            for j in range(n_gauss_psf_i):
                weight_ij, sigma_ij = psf_i['weight'][j], psf_i['sigma'][j]
                label = f'[weight, sigma of comp {j+1} of psf {i+1}]'
                line = f"{weight_ij} {sigma_ij}{tab}{label}\n"
                f.write(line)
            f.write(f'"infil/{kin_i.binfile}"\n')
        f.write('"datfil/mass_aper.dat"')
        f.close()

    def get_orbit_ics(self):
        """Execute the bash script to calculate orbit ICs
        """
        cur_dir = os.getcwd()
        os.chdir(self.mod_dir)
        cmdstr = self.write_executable_for_ics()
        self.logger.info('Calculating initial conditions')
        # p = subprocess.call('bash '+cmdstr, shell=True)
        p = subprocess.run('bash '+cmdstr,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,
                           shell=True)
        log_file = f'Logfile: {self.mod_dir}datfil/orbstart.log.'
        if not p.stdout.decode("UTF-8"):
            self.logger.info(f'...done - {cmdstr} exit code {p.returncode}. '
                             f'{log_file}')
        else:
            text = f'{cmdstr} exit code {p.returncode}. ERROR. ' \
                   f'Message: {p.stdout.decode("UTF-8")}{log_file}'
            self.logger.error(text)
            raise RuntimeError(text)
        os.chdir(cur_dir)

    def write_executable_for_ics(self):
        """Write the bash script to calculate orbit ICs
        """
        cmdstr = 'cmd_orb_start'
        #create the fortran executable
        txt_file = open(cmdstr, "w")
        txt_file.write('#!/bin/bash' + '\n')
        tmp = '/orbitstart < infil/orbstart.in >> datfil/orbstart.log\n'
        txt_file.write(f'{self.legacy_directory}{tmp}')
        txt_file.close()
        # the name of the executable must be returned to use in subprocess.call
        return cmdstr

    def get_orbit_library_par(self):
        """Execute the bash script to calculate orbit libraries in parallel
        """
        # move to model directory
        cur_dir = os.getcwd()
        os.chdir(self.mod_dir)
        cmdstr = self.write_executable_for_integrate_orbits_par()
        self.logger.info('Integrating orbit library tube and box orbits')
        # p = subprocess.call('bash '+cmdstr_tube, shell=True)
        p = subprocess.run('bash '+cmdstr,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,
                           shell=True)
        log_files = f'Logfiles: {self.mod_dir}datfil/orblib.log, ' \
                    f'{self.mod_dir}datfil/orblibbox.log, ' \
                    f'{self.mod_dir}datfil/triaxmass.log, ' \
                    f'{self.mod_dir}datfil/triaxmassbin.log.'
        if not p.stdout.decode("UTF-8"):
            self.logger.info(f'...done - {cmdstr} exit code '
                             f'{p.returncode}. {log_files}')
        else:
            text=f'{cmdstr} exit code {p.returncode}. ERROR. ' \
                 f'Message: {p.stdout.decode("UTF-8")}{log_files}'
            self.logger.error(text)
            raise RuntimeError(text)
        # move back to original directory
        os.chdir(cur_dir)

    def get_orbit_library(self):
        """Execute the bash script to calculate orbit libraries
        """
        # move to model directory
        cur_dir = os.getcwd()
        os.chdir(self.mod_dir)
        cmdstr_tube, cmdstr_box = self.write_executable_for_integrate_orbits()
        self.logger.info('Integrating orbit library tube orbits')
        # p = subprocess.call('bash '+cmdstr_tube, shell=True)
        p = subprocess.run('bash '+cmdstr_tube,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,
                           shell=True)
        log_files = f'Logfiles: {self.mod_dir}datfil/orblib.log, ' \
                    f'{self.mod_dir}datfil/triaxmass.log, ' \
                    f'{self.mod_dir}datfil/triaxmassbin.log.'
        if not p.stdout.decode("UTF-8"):
            self.logger.info(f'...done - {cmdstr_tube} exit code '
                             f'{p.returncode}. {log_files}')
        else:
            text=f'{cmdstr_tube} exit code {p.returncode}. ERROR. ' \
                 f'Message: {p.stdout.decode("UTF-8")}{log_files}'
            self.logger.error(text)
            raise RuntimeError(text)
        self.logger.info('Integrating orbit library box orbits')
        # p = subprocess.call('bash '+cmdstr_box, shell=True)
        p = subprocess.run('bash '+cmdstr_box,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,
                           shell=True)
        log_file = f'Logfile: {self.mod_dir}datfil/orblibbox.log.'
        if not p.stdout.decode("UTF-8"):
            self.logger.info(f'...done - {cmdstr_box} exit code '
                             f'{p.returncode}. {log_file}')
        else:
            text = f'{cmdstr_box} exit code {p.returncode}. ERROR. ' \
                   f'Message: {p.stdout.decode("UTF-8")}{log_file}'
            self.logger.error(text)
            raise RuntimeError(text)
        # move back to original directory
        os.chdir(cur_dir)

    def write_executable_for_integrate_orbits_par(self):
        """Write the bash script to calculate orbit libraries
        """
        if self.settings['use_new_mirroring']:
            orb_prgrm = 'orblib_new_mirror'
        else:
            orb_prgrm = 'orblib'
        cmd_string = 'cmd_tube_box_orbs'
        txt_file = open(cmd_string, "w")
        txt_file.write('#!/bin/bash\n')
        txt_file.write('(rm -f datfil/orblib.dat.tmp datfil/orblib.dat\n')
        txt_file.write(f'{self.legacy_directory}/{orb_prgrm} < infil/orblib.in '
                        '>> datfil/orblib.log\n')
        txt_file.write('rm -f datfil/mass_qgrid.dat datfil/mass_radmass.dat '
                        'datfil/mass_aper.dat\n')
        txt_file.write(f'{self.legacy_directory}/triaxmass '
                        '< infil/triaxmass.in >> datfil/triaxmass.log\n')
        txt_file.write(f'{self.legacy_directory}/triaxmassbin '
                        '< infil/triaxmassbin.in >> datfil/triaxmassbin.log\n')
        txt_file.write('rm -f datfil/orblib.dat.bz2 '
                        '&& bzip2 -k datfil/orblib.dat\n')
        txt_file.write('rm datfil/orblib.dat) &\n')
        txt_file.write('orblib=$!\n')
        txt_file.write('(rm -f datfil/orblibbox.dat.tmp datfil/orblibbox.dat\n')
        txt_file.write(f'{self.legacy_directory}/{orb_prgrm} '
                       '< infil/orblibbox.in >> datfil/orblibbox.log\n')
        txt_file.write('rm -f datfil/orblibbox.dat.bz2 '
                       '&& bzip2 -k datfil/orblibbox.dat\n')
        txt_file.write('rm datfil/orblibbox.dat) &\n')
        txt_file.write('orblibbox=$!\n')
        txt_file.write('wait $orblib $orblibbox\n')
        txt_file.close()
        # returns the name of the executables
        return cmd_string

    def write_executable_for_integrate_orbits(self):
        """Write the bash script to calculate orbit libraries
        """
        if self.settings['use_new_mirroring']:
            orb_prgrm = 'orblib_new_mirror'
        else:
            orb_prgrm = 'orblib'
        # tubeorbits
        cmdstr_tube = 'cmd_tube_orbs'
        txt_file = open(cmdstr_tube, "w")
        txt_file.write('#!/bin/bash\n')
        txt_file.write('rm -f datfil/orblib.dat.tmp datfil/orblib.dat '
                       'datfil/orblib.dat.bz2\n')
        txt_file.write(f'{self.legacy_directory}/{orb_prgrm} < infil/orblib.in '
                       '>> datfil/orblib.log\n')
        txt_file.write('rm -f datfil/mass_qgrid.dat datfil/mass_radmass.dat '
                       'datfil/mass_aper.dat\n')
        txt_file.write(f'{self.legacy_directory}/triaxmass '
                       '< infil/triaxmass.in >> datfil/triaxmass.log\n')
        txt_file.write(f'{self.legacy_directory}/triaxmassbin '
                       '< infil/triaxmassbin.in >> datfil/triaxmassbin.log\n')
        txt_file.write(
            'bzip2 -kc datfil/orblib.dat > datfil/orblib.dat.staging.bz2 '
            '&& mv datfil/orblib.dat.staging.bz2 datfil/orblib.dat.bz2\n')
        txt_file.write('rm datfil/orblib.dat\n')
        txt_file.close()
        # boxorbits
        cmdstr_box = 'cmd_box_orbs'
        txt_file = open(cmdstr_box, "w")
        txt_file.write('#!/bin/bash\n')
        txt_file.write('rm -f datfil/orblibbox.dat.tmp datfil/orblibbox.dat '
                       'datfil/orblibbox.dat.bz2\n')
        txt_file.write(f'{self.legacy_directory}/{orb_prgrm} '
                       '< infil/orblibbox.in >> datfil/orblibbox.log\n')
        txt_file.write(
            'bzip2 -kc datfil/orblibbox.dat > datfil/orblibbox.dat.staging.bz2 '
            '&& mv datfil/orblibbox.dat.staging.bz2 datfil/orblibbox.dat.bz2\n')
        txt_file.write('rm datfil/orblibbox.dat\n')
        txt_file.close()
        # returns the name of the executables
        return cmdstr_tube, cmdstr_box

    def read_ics(self):
        # ...
        pass

    def read_orbit_base(self, fileroot):
        """
        Read orbit library from file datfil/{fileroot}.dat.bz2'

        Parameters
        ----------
        fileroot : string
            this will probably be either 'orblib' or 'orblibbox'

        Returns
        -------
        Histogram
            the orbit library stored in a Histogram object

        """
        cur_dir = os.getcwd()
        os.chdir(self.mod_dir)
        # unzip orblib to a temproary file with ml value attached
        # ml value is needed to prevent different processes clashing
        ml = self.parset['ml']
        aaa = subprocess.run(['bunzip2', '-c', f'datfil/{fileroot}.dat.bz2'],
                             stdout=subprocess.PIPE)
        tmpfname = f'datfil/{fileroot}_{ml}.dat'
        tmpf = open(tmpfname, "wb")
        tmpf.write(aaa.stdout)
        tmpf.close()
        # read the fortran file
        orblibf = FortranFile(tmpfname, 'r')
        # remove temproary file
        subprocess.call(['rm', tmpfname])
        # read size of orbit library
        # from integrator_setup_write, lines 506 - 5129:
        tmp = orblibf.read_ints(np.int32)
        norb, t2, t3, t4, ndith = tmp
        # from qgrid_setup_write, lines 2339-1350:
        quad_light_grid_sizes = orblibf.read_ints(np.int32)
        size_ql, size_qph, size_qth, size_qlr = quad_light_grid_sizes
        quad_lr = orblibf.read_reals(float)
        quad_lth = orblibf.read_reals(float)
        quad_lph = orblibf.read_reals(float)
        # from histogram_setup_write, lines 1917-1926:
        tmp = orblibf.read_record(np.int32, np.int32, float)
        nconstr = tmp[0][0] # = total number of apertures for ALL kinematics
        nvhist = tmp[1][0] # = (nvbins-1)/2 for histo of FIRST kinematic set
        dvhist = tmp[2][0] # = delta_v in histogram for FIRST kinematic set
        # these nvhist and dvhist are for the first kinematic set only
        # however, orbit are stored N times where N = number of kinematic sets
        # histogram settings for other N-1 sets may be different from the first
        # these aren't stored in orblib.dat so must read from kinematics objects
        stars = self.system.get_component_from_class(
            physys.TriaxialVisibleComponent
            )
        n_kins = len(stars.kinematic_data)
        hist_widths = [k.hist_width for k in stars.kinematic_data]
        hist_centers = [k.hist_center for k in stars.kinematic_data]
        hist_bins = [k.hist_bins for k in stars.kinematic_data]
        self.logger.debug('Checking number of velocity bins...')
        error_msg = 'must have odd number of velocity bins for all kinematics'
        assert np.all(np.array(hist_bins) % 1==0), error_msg
        self.logger.debug('...checks ok.')
        n_apertures = [len(k.data) for k in stars.kinematic_data]
        # get index linking  kinematic set to aperture
        # kin_idx_per_ap[i] = N <--> aperture i is from kinematic set N
        kin_idx_per_ap = [np.zeros(n_apertures[i], dtype=int)+i
                          for i in range(n_kins)]
        kin_idx_per_ap = np.concatenate(kin_idx_per_ap)
        kin_idx_per_ap = np.array(kin_idx_per_ap, dtype=int)
        # below we loop i_ap from 1-n_total_apertures but will need the index of
        # i_ap for the relevant kinematic set: we use `idx_ap_reset` to do this
        cum_n_apertures = np.cumsum(n_apertures)
        idx_ap_reset = np.concatenate(([0], cum_n_apertures[:-1]))
        # set up a list of arrays to hold the results
        tmp = zip(hist_bins,n_apertures)
        velhist0 = [np.zeros((norb, nv, na)) for (nv,na) in tmp]
        # Next read the histograms themselves.
        orbtypes = np.zeros((norb, ndith**3), dtype=int)
        nbins_vhist = 2*nvhist + 1
        velhist = np.zeros((norb, nbins_vhist, nconstr))
        density_3D = np.zeros((norb, size_qlr, size_qth, size_qph))
        for j in range(norb):
            t1,t2,t3,t4,t5 = orblibf.read_ints(np.int32)
            orbtypes[j, :] = orblibf.read_ints(np.int32)
            quad_light = orblibf.read_reals(float)
            quad_light = np.reshape(quad_light, quad_light_grid_sizes[::-1])
            # quad_light stores orbit features in 3D (r,th,phi) bins.
            # Quad_light[ir,it,ip,XXX] stores
            # - the zeroth moment i.e. density for XXX=0,
            # - the first moments x,y,z,vx,vy,vz for XXX=slice(1,7)
            # - 2nd moments vx^2,vy^2,vz^2,vx*vy,vy*vz,vz*vx for XXX=slice(7,13)
            # - an averaged orbit classification for XXX=slice(13,None)
            # in the bin indexed by (ir,it,ip).
            # We need to extract 3D density for use in weight solving.
            density_3D[j] = quad_light[:,:,:,0]
            for i_ap in range(nconstr):
                kin_idx = kin_idx_per_ap[i_ap]
                i_ap0 = i_ap - idx_ap_reset[kin_idx]
                ivmin, ivmax = orblibf.read_ints(np.int32)
                if ivmin <= ivmax:
                    nv0 = (hist_bins[kin_idx]-1)/2
                    # ^--- this is an integer since hist_bins is odd
                    nv0 = int(nv0)
                    tmp = orblibf.read_reals(float)
                    velhist0[kin_idx][j, ivmin+nv0:ivmax+nv0+1, i_ap0] = tmp
        orblibf.close()
        os.chdir(cur_dir)
        velhists = []
        for i in range(n_kins):
            center0 = hist_centers[i]
            width0 = hist_widths[i]
            bins0 = hist_bins[i]
            idx_center = (bins0-1)/2 # this is an integer since hist_bins is odd
            idx_center = int(idx_center)
            dvhist0 = width0/bins0
            vedg = np.arange(bins0+1) * dvhist0
            v = (vedg[1:] + vedg[:-1])/2.
            v_cent = v[idx_center]
            delta_v = center0 - v_cent
            vedg -= v_cent
            vvv = dyn_kin.Histogram(xedg=vedg,
                                    y=velhist0[i],
                                    normalise=False)
            velhists += [vvv]
        return velhists, density_3D

    def duplicate_flip_and_interlace_orblib(self, orblib):
        """mirror the tube orbits

        Take an orbit library, create a duplicate library with the velocity
        signs flipped, then interlace the two i.e. so that resulting library
        alternates between flipped/unflipped. This creates an orbit library
        consistent with the Fortran output, enforcing the ordering created by
        the for loops in lines 157-178 of triaxnnls_CRcut.f90

        Parameters
        ----------
        orblib : ``dyn.kinematics.Histogram``

        Returns
        -------
        ``dyn.kinematics.Histogram``
            the duplicated, flipped and interlaced orblib

        """
        self.logger.debug('Checking for symmetric velocity array...')
        error_msg = 'velocity array must be symmetric'
        assert np.allclose(orblib.xedg, -orblib.xedg[::-1]), error_msg
        self.logger.debug('...check ok.')
        losvd = orblib.y
        n_orbs, n_vel_bins, n_spatial_bins = losvd.shape
        reveresed_losvd = losvd[:, ::-1, :]
        new_losvd = np.zeros((2*n_orbs, n_vel_bins, n_spatial_bins))
        new_losvd[0::2] = losvd
        new_losvd[1::2, :] = reveresed_losvd
        new_orblib = dyn_kin.Histogram(xedg=orblib.xedg,
                                       y=new_losvd,
                                       normalise=False)
        return new_orblib

    def combine_orblibs(self, orblib1, orblib2):
        """Combine two LOSVD histograms into one.

        Parameters
        ----------
        orblib1 : ``dyn.kinematics.Histogram``
        orblib2 : ``dyn.kinematics.Histogram``

        Returns
        -------
        ``dyn.kinematics.Histogram``
            the combined orbit libraries

        """
        # check orblibs are compatible
        n_orbs1, n_vel_bins1, n_spatial_bins1 = orblib1.y.shape
        n_orbs2, n_vel_bins2, n_spatial_bins2 = orblib2.y.shape
        self.logger.debug('Checking number of velocity bins...')
        error_msg = 'orblibs have different number of velocity bins'
        assert n_vel_bins1==n_vel_bins2, error_msg
        self.logger.debug('Checking velocity arrays...')
        error_msg = 'orblibs have different velocity arrays'
        assert np.array_equal(orblib1.x, orblib2.x), error_msg
        self.logger.debug('Checking number of spatial bins...')
        error_msg = 'orblibs have different number of spatial bins'
        assert n_spatial_bins1==n_spatial_bins2, error_msg
        self.logger.debug('...checks ok.')
        new_losvd = np.zeros((n_orbs1 + n_orbs2,
                              n_vel_bins1,
                              n_spatial_bins1))
        new_losvd[:n_orbs1] = orblib1.y
        new_losvd[n_orbs1:] = orblib2.y
        new_orblib = dyn_kin.Histogram(xedg=orblib1.xedg,
                                       y=new_losvd,
                                       normalise=False)
        return new_orblib

    def read_losvd_histograms(self):
        """Read the orbit library

        Read box orbits and tube orbits, mirrors the latter, and combines.
        Rescales the velocity axis according to the ``ml`` value. Sets LOSVDs
        and 3D grid/aperture masses of the combined orbit library.

        Returns
        -------
        Sets the attributes:
            -   ``self.losvd_histograms``: a list, whose i'th entry is a
                ``dyn.kinematics.Histogram`` object holding the orbit lib LOSVDs
                binned for the i'th kinematic set
            -   ``self.intrinsic_masses``: 3D grid/intrinsic masses of orbit lib
            -   ``self.projected_masses``: aperture/proj. masses of orbit lib
            -   ``self.n_orbs``: number of orbits in the orbit library

        """
        # TODO: check if this ordering is compatible with weights read in by
        # LegacyWeightSolver.read_weights
        tube_orblib, tube_density_3D = self.read_orbit_base('orblib')
        # tube orbits are mirrored/flipped and used twice
        tmp = []
        for tube_orblib0 in tube_orblib:
            tmp += [self.duplicate_flip_and_interlace_orblib(tube_orblib0)]
        tube_orblib = tmp
        tube_density_3D = np.repeat(tube_density_3D, 2, axis=0)
        # read box orbits
        box_orblib, box_density_3D = self.read_orbit_base('orblibbox')
        # combine orblibs
        orblib = []
        for (t0, b0) in zip(tube_orblib, box_orblib):
            orblib0 = self.combine_orblibs(t0, b0)
            orblib += [orblib0]
        # combine density_3D arrays
        density_3D = np.vstack((tube_density_3D, box_density_3D))
        ml_current = self.parset['ml']
        ml_original = self.get_ml_of_original_orblib()
        scale_factor = np.sqrt(ml_current/ml_original)
        nkins = len(orblib)
        for i in range(nkins):
            orblib[i].scale_x_values(scale_factor)
        self.losvd_histograms = orblib
        self.intrinsic_masses = density_3D
        self.n_orbs = self.losvd_histograms[0].y.shape[0]
        proj_mass = [np.sum(self.losvd_histograms[i].y,1) for i in range(nkins)]
        self.projected_masses = proj_mass

    def get_ml_of_original_orblib(self):
        """Get ``ml`` of original orblib with shared parameters

        The original ``ml`` is required to rescale orbit libraries for rescaled
        potentials. This method reads it from the first entry of 9th from bottom
        line of the file ``infil/parameters_pot.in``

        Returns
        -------
        float
            the original ``ml``

        """
        infile = self.mod_dir + 'infil/parameters_pot.in'
        lines = [line.rstrip('\n').split() for line in open(infile)]
        ml_original = float((lines[-9])[0])
        return ml_original

    def read_orbit_property_file_base(self, file, ncol, nrow):
        """
        read in `datfil/orblib.dat_orbclass.out`

        which stores the information of all the orbits stored in the orbit library
        norb = nrow = nE * nI2 * nI3 * ndithing^3
        for each orbit, the time averaged values are stored:
        lx, ly ,lz, r = sum(sqrt( average(r^2) )), Vrms^2 = average(vx^2 + vy^2 + vz^2 + 2vx*vy + 2vxvz + 2vxvy)

        The file was stored by the fortran code orblib_f.f90 integrator_find_orbtype
        """
        data=[]
        lines = [line.rstrip('\n').split() for line in open(file)]
        i = 0
        while i < len(lines):
            for x in lines[i]:
                data.append(np.double(x))
            i += 1
        data=np.array(data)
        data=data.reshape((5,ncol,nrow), order='F')
        return data

    def read_orbit_property_file(self):
        """Read the file `*orbclass.out` files

        This file contains angular momenta Lx, Ly, Lz and v_rms of the orbits.
        Results are stored in and `self.table`

        """
        nE = self.settings['nE']
        nI2 = self.settings['nI2']
        nI3 = self.settings['nI3']
        ndither = self.settings['dithering']
        norb = int(nE*nI2*nI3)
        ncol = ndither**3
        nrow = norb
        orbclass1 = self.read_orbit_property_file_base(
            self.mod_dir+'datfil/orblib.dat_orbclass.out',
            ncol=ncol,
            nrow=nrow)
        orbclass2 = self.read_orbit_property_file_base(
            self.mod_dir+'datfil/orblibbox.dat_orbclass.out',
            ncol=ncol,
            nrow=nrow)
        orbclass=np.dstack((orbclass1,orbclass1,orbclass2))
        orbclass1a=np.copy(orbclass1)
        orbclass1a[0:3,:,:] *= -1
        for i in range(int(0), norb):
            orbclass[:,:,i*2]=orbclass1[:, :, i]
            orbclass[:,:,i*2 + 1]=orbclass1a[:, :, i]
        orb_properties = table.QTable()
        orb_properties['Lx'] = orbclass[0,:,:].T * u.km * u.km/u.s
        orb_properties['Ly'] = orbclass[1,:,:].T * u.km * u.km/u.s
        orb_properties['Lz'] = orbclass[2,:,:].T * u.km * u.km/u.s
        orb_properties['r'] = orbclass[3,:,:].T * u.km
        orb_properties['Vrms'] = orbclass[4,:,:].T**0.5 * u.km/u.s
        r_vrms = orb_properties['r']*orb_properties['Vrms']
        orb_properties['lmd_x'] = orb_properties['Lx']/r_vrms
        orb_properties['lmd_y'] = orb_properties['Ly']/r_vrms
        orb_properties['lmd_z'] = orb_properties['Lz']/r_vrms
        orb_properties['r'] = orb_properties['r'].to(u.kpc)
        orb_properties['L'] = (
            orb_properties['Lx']**2 + 
            orb_properties['Ly']**2 + 
            orb_properties['Lz']**2)**0.5
        orb_properties['lmd'] = (orb_properties['L']/r_vrms).to(u.dimensionless_unscaled)
        self.orb_properties = orb_properties

    def find_threshold_angular_momentum(self, make_diagnostic_plots=False):
        orb_properties = self.orb_properties
        dl = np.linspace(-1, 1, 12)[6]
        dl *= 1.
        if make_diagnostic_plots:
            kw_hist = {'range':(-1,1), 'bins':11, 'alpha':0.3}
            _ = plt.hist(np.ravel(orb_properties['Lx'].value/1e18), **kw_hist)
            _ = plt.hist(np.ravel(orb_properties['Ly'].value/1e18), **kw_hist)
            _ = plt.hist(np.ravel(orb_properties['Lz'].value/1e18), **kw_hist)
            plt.axvline(dl, ls=':', color='k')
            plt.axvline(-dl, ls=':', color='k')
        return dl

    def classify_orbits(self, make_diagnostic_plots=False):
        orb_properties = self.orb_properties
        dl = self.find_threshold_angular_momentum(
            make_diagnostic_plots=make_diagnostic_plots
            )
        # find box orbits
        bool_box = (
            (np.abs(orb_properties['Lx'].value/1e18) <= dl) &
            (np.abs(orb_properties['Ly'].value/1e18) <= dl) &
            (np.abs(orb_properties['Lz'].value/1e18) <= dl)
        )
        idx_box = np.where(bool_box)
        # find "true" tube orbits i.e. with exactly one component of L =/= 0
        bool_xtube = (
            (np.abs(orb_properties['Lx'].value/1e18) > dl) &
            (np.abs(orb_properties['Ly'].value/1e18) <= dl) &
            (np.abs(orb_properties['Lz'].value/1e18) <= dl)
        )
        idx_xtube = np.where(bool_xtube)
        bool_ytube = (
            (np.abs(orb_properties['Lx'].value/1e18) <= dl) &
            (np.abs(orb_properties['Ly'].value/1e18) > dl) &
            (np.abs(orb_properties['Lz'].value/1e18) <= dl)
        )
        idx_ytube = np.where(bool_ytube)
        bool_ztube = (
            (np.abs(orb_properties['Lx'].value/1e18) <= dl) &
            (np.abs(orb_properties['Ly'].value/1e18) <= dl) &
            (np.abs(orb_properties['Lz'].value/1e18) > dl)
        )
        idx_ztube = np.where(bool_ztube)
        # find tube-ish orbits i.e. with one component of L larger than other 2
        bool_xtish = (
            (bool_box==False)&
            (np.abs(orb_properties['Lx']) > np.abs(orb_properties['Ly'])) &
            (np.abs(orb_properties['Lx']) > np.abs(orb_properties['Lz']))
        )
        bool_ytish = (
            (bool_box==False) &
            (np.abs(orb_properties['Ly']) > np.abs(orb_properties['Lx'])) &
            (np.abs(orb_properties['Ly']) > np.abs(orb_properties['Lz']))
        )
        bool_ztish = (
            (bool_box==False) &
            (np.abs(orb_properties['Lz']) > np.abs(orb_properties['Lx'])) &
            (np.abs(orb_properties['Lz']) > np.abs(orb_properties['Ly']))
        )
        # find any remaining orbits
        bool_other = (
            (bool_box==False) &
            (bool_xtish==False) &
            (bool_ytish==False) &
            (bool_ztish==False))
        # log
        n_orb_tot = bool_box.size
        n_box = np.sum(bool_box)
        n_xtish = np.sum(bool_xtish)
        n_ytish = np.sum(bool_ytish)
        n_ztish = np.sum(bool_ztish)
        n_other = np.sum(bool_other)
        n_sum = n_box + n_xtish + n_ytish + n_ztish + n_other
        n_xt_exact = np.sum(bool_xtube)
        n_yt_exact = np.sum(bool_ytube)
        n_zt_exact = np.sum(bool_ztube)
        def percent(x):
            return f'{100.*x:.1f}%'
        self.logger.info('Orbit library classification:')
        self.logger.info(f'    - {percent(n_box/n_orb_tot)} box')
        self.logger.info(f'    - {percent(n_xtish/n_orb_tot)} x-tubes')
        self.logger.info(f'    - {percent(n_ytish/n_orb_tot)} y-tubes')
        self.logger.info(f'    - {percent(n_ztish/n_orb_tot)} z-tubes')
        self.logger.info(f'    - {percent(n_other/n_orb_tot)} other types')
        self.logger.info('Amongst tubes, % with only one nonzero component of L:')
        self.logger.info(f'    - {percent(n_xt_exact/n_xtish)} of x-tubes')
        self.logger.info(f'    - {percent(n_yt_exact/n_ytish)} of y-tubes')
        self.logger.info(f'    - {percent(n_zt_exact/n_ztish)} of z-tubes')
        self.logger.info('Orbit library classification DONE.')
        # save the output
        orb_classification = {
            'bool_box':bool_box,
            'bool_xtish':bool_xtish,
            'bool_ytish':bool_ytish,
            'bool_ztish':bool_ztish,
            'bool_other':bool_other
        }
        self.orb_classification = orb_classification

    def get_projection_tensor(self, minr=None, maxr=None, nr=50, nl=61):
        # skip if this method has been run before with the same input
        projection_tensor_pars = {'minr':minr, 
                                  'maxr':maxr, 
                                  'nr':nr,
                                  'nl':nl}
        if hasattr(self, 'projection_tensor_pars'):
            if self.projection_tensor_pars == projection_tensor_pars:
                return
            else:
                self.projection_tensor_pars = projection_tensor_pars
        # otherwise, continue...
        if hasattr(self, 'orb_properties') == False:
            self.read_orbit_property_file()
        orb_properties = self.orb_properties
        if hasattr(self, 'orb_classification') == False:
            self.classify_orbits()
        if minr is None:
            minr = np.min(orb_properties['r']).value
        if maxr is None:
            maxr = np.max(orb_properties['r']).value
        log10_r_rng = (np.log10(minr), np.log10(maxr))
        lmd_rng = (-1, 1)
        tot_lmd_rng = (0, 1)
        # store ranges for use in plots
        self.projection_tensor_rng = {
            'log10_r_rng':log10_r_rng,
            'lmd_rng':lmd_rng,
            'tot_lmd_rng':tot_lmd_rng,
        }
        # store ranges for use in plots
        log10_r_edg = np.linspace(*log10_r_rng, nr+1)
        lmd_edg = np.linspace(*lmd_rng, nl+1)
        tot_lmd_edg = np.linspace(*tot_lmd_rng, nl+1)
        r_idx = np.digitize(np.log10(orb_properties['r'].value), bins=log10_r_edg)
        lmd_x_idx = np.digitize(orb_properties['lmd_x'].value, bins=lmd_edg)
        lmd_y_idx = np.digitize(orb_properties['lmd_y'].value, bins=lmd_edg)
        lmd_z_idx = np.digitize(orb_properties['lmd_z'].value, bins=lmd_edg)
        tot_lmd_idx = np.digitize(orb_properties['lmd'].value, bins=tot_lmd_edg)
        bundle_idx, orbit_idx = np.indices(r_idx.shape)
        # make (sparse matrix representation of) projection tensor
        projection = []
        orbtypes = ['xtish', 'ytish', 'ztish', 'box']
        ang_mom_indixes = [lmd_x_idx, lmd_y_idx, lmd_z_idx, tot_lmd_idx]
        for (orbtype00, l_idx00) in zip(orbtypes, ang_mom_indixes):
            str00 = f'bool_{orbtype00}'
            bool00 = self.orb_classification[str00]
            # decrease r_idx/L_idx by 1 so they are 0-index
            coords = np.array([bundle_idx[bool00], 
                               orbit_idx[bool00], 
                               r_idx[bool00]-1, 
                               l_idx00[bool00]-1])
            # remove entries where orbit is outside bounds
            idx_keep = np.where(
                (coords[2,:]>-1) & (coords[3,:]>-1) & (coords[2,:]<nr) & (coords[3,:]<nl)
                )
            coords = coords[:, idx_keep[0]]
            # create binary sparse matrix, with 1 entry per particle per bundle
            projection00 = sparse.COO(coords, 1, shape=r_idx.shape+(nr,nl))
            # average over individual orbits in a bundle
            projection00 = np.mean(projection00, 1)
            projection += [projection00]
        projection = np.stack(projection)
        projection = np.moveaxis(projection, 1, 3)
        self.projection_tensor = projection

# end
