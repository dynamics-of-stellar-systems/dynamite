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
from dynamite.constants import PARSEC_KM

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
        if len(config.all_models.table) == 0:
            self.velocity_scaling_factor = 1.0
        else:
            mod = config.all_models.get_model_from_parset(self.parset)
            self.velocity_scaling_factor = \
                config.all_models.get_model_velocity_scaling_factor(model=mod)

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
            if self.system.is_bar_disk_system():
                stars = self.system.get_unique_bar_component()
            else:
                stars = self.system.get_unique_triaxial_visible_component()
            # create the kinematics and populations input files for each
            # kinematic dataset and population dataset with own apertures
            kin_pops = stars.kinematic_data
            kin_pops += [p for p in stars.population_data if p.kin_aper is None]
            for data_set in kin_pops:
                # copy aperture and bins files across
                shutil.copyfile(self.in_dir + data_set.aperturefile,
                                self.mod_dir + f'infil/{data_set.aperturefile}')
                shutil.copyfile(self.in_dir + data_set.binfile,
                                self.mod_dir + f'infil/{data_set.binfile}')
            # calculate orbit libary
            file1 = 'begin.dat'
            file2 = 'beginbox.dat'
            check1 = os.path.isfile(self.mod_dir + f'datfil/{file1}')
            check2 = os.path.isfile(self.mod_dir + f'datfil/{file2}')
            if check1 + check2 != 2:
                if check1:
                    os.remove(self.mod_dir + f'datfil/{file1}')
                if check2:
                    os.remove(self.mod_dir + f'datfil/{file2}')
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
        if self.system.is_bar_disk_system():
            stars = self.system.get_unique_bar_component()
        else:
            stars = self.system.get_unique_triaxial_visible_component()
        bh = self.system.get_component_from_class(physys.Plummer)
        # used to derive the viewing angles
        if self.system.is_bar_disk_system():
            if self.system.is_bar_disk_system_with_angles():
                theta = self.parset[f'theta-{stars.name}']
                phi = self.parset[f'phi-{stars.name}']
                psi = self.parset[f'psi-{stars.name}']
            else:
                q = self.parset[f'q-{stars.name}']
                p = self.parset[f'p-{stars.name}']
                u = self.parset[f'u-{stars.name}']
                qdisk = self.parset[f'qdisk-{stars.name}']
                theta, psi, phi = stars.triax_pqu2tpp_bar(p,q,u,qdisk)
                phi = -phi ## FIX ME
        else:
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
        if self.system.is_bar_disk_system():
            len_disk_pot = len(stars.disk_pot.data)
            header_string_pot = str(len_mge_pot + len_disk_pot) + " 1 " + str(len_mge_pot) + " " + str(len_disk_pot)
            len_disk_lum = len(stars.disk_lum.data)
            header_string_lum = str(len_mge_lum + len_disk_lum) + " 1 " + str(len_mge_lum) + " " + str(len_disk_lum)
            text += f"\n{self.parset['omega']}"
            mge_pot = stars.mge_pot + stars.disk_pot
            mge_lum = stars.mge_lum + stars.disk_lum
        else:
            header_string_pot = str(len_mge_pot)
            header_string_lum = str(len_mge_lum)
            mge_pot = stars.mge_pot
            mge_lum = stars.mge_lum

        # parameters_pot.in
        np.savetxt(path + 'parameters_pot.in',
                   mge_pot.data,
                   header=header_string_pot,
                   footer=text,
                   comments='',
                   fmt=['%10.2f','%10.5f','%10.5f','%10.2f'])
        # parameters_lum.in
        np.savetxt(path + 'parameters_lum.in',
                   mge_lum.data,
                   header=header_string_lum,
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
        n_psf_kin = len(stars.kinematic_data)
        psf_pop_idx = [i for i, p in enumerate(stars.population_data)
                       if p.kin_aper is None]  # pops with their own apertures
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
            line = f"{self.settings['orbital_periods']}{tab}{label}\n"
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
            label = '[number of psfs of the kinematic data + population data]'
            line = f"{n_psf_kin+len(psf_pop_idx)}{tab}{label}\n"
            f.write(line)
            for i in range(n_psf_kin):
                psf_i = stars.kinematic_data[i].PSF
                n_gauss_psf_i = len(psf_i['sigma'])
                label = f'[# of gaussians in kinematics psf {i+1}]'
                line = f"{n_gauss_psf_i}{tab}{label}\n"
                f.write(line)
            for i in psf_pop_idx:
                psf_i = stars.population_data[i].PSF
                n_gauss_psf_i = len(psf_i['sigma'])
                label = f'[# of gaussians in populations psf {i+1}]'
                line = f"{n_gauss_psf_i}{tab}{label}\n"
                f.write(line)
            for i in range(n_psf_kin):
                psf_i = stars.kinematic_data[i].PSF
                for j in range(len(psf_i['sigma'])):
                    weight_ij, sigma_ij = psf_i['weight'][j], psf_i['sigma'][j]
                    label = f'[weight, sigma of comp {j+1} of kins psf {i+1}]'
                    line = f"{weight_ij} {sigma_ij}{tab[:-1]}{label}\n"
                    f.write(line)
            for i in psf_pop_idx:
                psf_i = stars.population_data[i].PSF
                for j in range(len(psf_i['sigma'])):
                    weight_ij, sigma_ij = psf_i['weight'][j], psf_i['sigma'][j]
                    label = f'[weight, sigma of comp {j+1} of pops psf {i+1}]'
                    line = f"{weight_ij} {sigma_ij}{tab[:-1]}{label}\n"
                    f.write(line)
            # every kinematic set and population set has a psf and aperture,
            # so n_psf_kin = n_aperture_kin and n_psf_pop = len(psf_pop_idx)
            # (note that some pops may re-use a kin psf, so not all pops
            # need to have a psf of their own)
            n_aperture_kin = n_psf_kin
            n_aperture_pop = len(psf_pop_idx)
            label = '[# of apertures of the kinematic data + population data]'
            line = f'{n_aperture_kin + n_aperture_pop}{tab}{label}\n'
            f.write(line)
            for i in np.arange(n_aperture_kin):
                kin_i = stars.kinematic_data[i]
                # aperturefile cant have a label since fortran reads whole line
                f.write(f'"infil/{kin_i.aperturefile}"\n')
                label = f'[use kinematics psf {i+1} for {kin_i.name}]'
                line = f'{i+1}{tab}{label}\n'
                f.write(line)
            for i0, i in enumerate(psf_pop_idx):  # i0 counts aperture sets
                pop_i = stars.population_data[i]
                # aperturefile cant have a label since fortran reads whole line
                f.write(f'"infil/{pop_i.aperturefile}"\n')
                label = f'[use populations psf {i+1} = ' \
                        f'total psf {i0+n_psf_kin+1} for {pop_i.name}]'
                line = f'{i0+n_psf_kin+1}{tab}{label}\n'
                f.write(line)
            for i in np.arange(n_aperture_kin):
                kin_i = stars.kinematic_data[i]
                label = '[vhist width, center and nbins for kinematics ' \
                        f'{kin_i.name}]'
                w, c, b = kin_i.hist_width, kin_i.hist_center, kin_i.hist_bins
                line = f'{w} {c} {b}{tab[:-2]}{label}\n'
                f.write(line)
            for i in psf_pop_idx:
                pop_i = stars.population_data[i]
                label = '[vhist width, center and nbins for populations ' \
                        f'{pop_i.name}]'
                w, c, b = pop_i.hist_width, pop_i.hist_center, pop_i.hist_bins
                line = f'{w} {c} {b}{tab[:-2]}{label}\n'
                f.write(line)
            for i in np.arange(n_aperture_kin):
                label = f'[use binning for kinematics aperture {1+i}? ' \
                        '0/1 = yes/no]'
                line = f'1{tab}{label}\n'
                f.write(line)
            for i in psf_pop_idx:
                label = f'[use binning for populations aperture {1+i}? ' \
                        '0/1 = yes/no]'
                line = f'1{tab}{label}\n'
                f.write(line)
            for i in np.arange(n_aperture_kin):
                kin_i = stars.kinematic_data[i]
                label = f'[binfile for kinematics aperture {1+i}]'
                line = f'"infil/{kin_i.binfile}"{tab[:-2]}{label}\n'
                f.write(line)
            for i in psf_pop_idx:
                pop_i = stars.population_data[i]
                label = f'[binfile for populations aperture {1+i}]'
                line = f'"infil/{pop_i.binfile}"{tab[:-3]}{label}\n'
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
        f = open(path + 'triaxmassbin.in', 'w')
        f.write('infil/parameters_lum.in\n')
        f.write(f'{n_psf_kin + len(psf_pop_idx)}{tab}'
                '[# of kinematics + populations apertures]\n')
        for i in range(n_psf_kin):
            kin_i = stars.kinematic_data[i]
            f.write(f'"infil/{kin_i.aperturefile}"\n')
            psf_i = kin_i.PSF
            n_gauss_psf_i = len(psf_i['sigma'])
            label = f'[# of gaussians in kinematics psf {i+1}]'
            line = f"{n_gauss_psf_i}{tab}{label}\n"
            f.write(line)
            for j in range(n_gauss_psf_i):
                weight_ij, sigma_ij = psf_i['weight'][j], psf_i['sigma'][j]
                label = f'[weight, sigma of comp {j+1} of kin psf {i+1}]'
                line = f"{weight_ij} {sigma_ij}{tab[:-1]}{label}\n"
                f.write(line)
            f.write(f'"infil/{kin_i.binfile}"\n')
        for i in psf_pop_idx:
            pop_i = stars.population_data[i]
            f.write(f'"infil/{pop_i.aperturefile}"\n')
            psf_i = pop_i.PSF
            n_gauss_psf_i = len(psf_i['sigma'])
            label = f'[# of gaussians in populations psf {i+1}]'
            line = f"{n_gauss_psf_i}{tab}{label}\n"
            f.write(line)
            for j in range(n_gauss_psf_i):
                weight_ij, sigma_ij = psf_i['weight'][j], psf_i['sigma'][j]
                label = f'[weight, sigma of comp {j+1} of pop psf {i+1}]'
                line = f"{weight_ij} {sigma_ij}{tab[:-1]}{label}\n"
                f.write(line)
            f.write(f'"infil/{pop_i.binfile}"\n')
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
            text = f'...failed! {cmdstr} exit code {p.returncode}. ' \
                   f'Message: {p.stdout.decode("UTF-8")}'
            if p.returncode == 127: # command not found
                text += 'Check DYNAMITE legacy_fortran executables.'
                self.logger.error(text)
                raise FileNotFoundError(text)
            else:
                text += f'{log_file} Be wary: DYNAMITE may crash...'
                self.logger.warning(text)
                raise RuntimeError(text)
        os.chdir(cur_dir)

    def write_executable_for_ics(self):
        """Write the bash script to calculate orbit ICs
        """
        cmdstr = 'cmd_orb_start'
        #create the fortran executable
        txt_file = open(cmdstr, "w")
        txt_file.write('#!/bin/bash' + '\n')
        if (self.system.is_bar_disk_system()):
            tmp = '/orbitstart_bar < infil/orbstart.in >> datfil/orbstart.log\n'
        else:
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
            text=f'...failed! {cmdstr} exit code {p.returncode}. ' \
                 f'Message: {p.stdout.decode("UTF-8")}'
            if p.returncode == 127: # command not found
                text += 'Check DYNAMITE legacy_fortran executables.'
                self.logger.error(text)
                raise FileNotFoundError(text)
            else:
                text += f'{log_files} Be wary: DYNAMITE may crash...'
                self.logger.warning(text)
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
            text=f'...failed! {cmdstr_tube} exit code {p.returncode}. ' \
                 f'Message: {p.stdout.decode("UTF-8")}'
            if p.returncode == 127: # command not found
                text += 'Check DYNAMITE legacy_fortran executables.'
                self.logger.error(text)
                raise FileNotFoundError(text)
            else:
                text += f'{log_files} Be wary: DYNAMITE may crash...'
                self.logger.warning(text)
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
            text = f'...failed! {cmdstr_box} exit code {p.returncode}. ' \
                   f'Message: {p.stdout.decode("UTF-8")}'
            if p.returncode == 127: # command not found
                text += 'Check DYNAMITE legacy_fortran executables.'
                self.logger.error(text)
                raise FileNotFoundError(text)
            else:
                text += f'{log_file} Be wary: DYNAMITE may crash...'
                self.logger.warning(text)
                raise RuntimeError(text)
        # move back to original directory
        os.chdir(cur_dir)

    def write_executable_for_integrate_orbits_par(self):
        """Write the bash script to calculate orbit libraries
        """
        if self.system.is_bar_disk_system():
            orb_prgrm = 'orblib_bar'
        elif self.settings['use_new_mirroring']:
            orb_prgrm = 'orblib_new_mirror'
        else:
            orb_prgrm = 'orblib'
        cmd_string = 'cmd_tube_box_orbs'
        txt_file = open(cmd_string, "w")
        txt_file.write('#!/bin/bash\n')
        txt_file.write('# first, check whether executables exist\n')
        for f_name in orb_prgrm, 'triaxmass', 'triaxmassbin':
            txt_file.write(f'test -e {self.legacy_directory}/{f_name} || ' +
                           f'{{ echo "File {self.legacy_directory}/{f_name} ' +
                           'not found." && exit 127; }\n')
        txt_file.write('(rm -f datfil/orblib.dat.tmp datfil/orblib.dat\n')
        txt_file.write(f'{self.legacy_directory}/{orb_prgrm} < infil/orblib.in '
                        '>> datfil/orblib.log\n')
        txt_file.write('rm -f datfil/mass_qgrid.dat datfil/mass_radmass.dat '
                        'datfil/mass_aper.dat\n')

        if self.system.is_bar_disk_system():
            txt_file.write(f'{self.legacy_directory}/triaxmass_bar '
                           '< infil/triaxmass.in >> datfil/triaxmass.log\n')
            txt_file.write(f'{self.legacy_directory}/triaxmassbin_bar '
                           '< infil/triaxmassbin.in >> datfil/triaxmassbin.log\n')
        else:
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
        if self.system.is_bar_disk_system():
            orb_prgrm = 'orblib_bar'
        elif self.settings['use_new_mirroring']:
            orb_prgrm = 'orblib_new_mirror'
        else:
            orb_prgrm = 'orblib'
        # tubeorbits
        cmdstr_tube = 'cmd_tube_orbs'
        txt_file = open(cmdstr_tube, "w")
        txt_file.write('#!/bin/bash\n')
        txt_file.write('# first, check whether executables exist\n')
        for f_name in orb_prgrm, 'triaxmass', 'triaxmassbin':
            txt_file.write(f'test -e {self.legacy_directory}/{f_name} || ' +
                           f'{{ echo "File {self.legacy_directory}/{f_name} ' +
                           'not found." && exit 127; }\n')
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
        txt_file.write('# first, check whether executable exists\n')
        txt_file.write(f'test -e {self.legacy_directory}/{orb_prgrm} || ' +
                       f'{{ echo "File {self.legacy_directory}/{orb_prgrm} ' +
                       'not found." && exit 127; }\n')
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

    def read_orbit_base(self, fileroot, return_instrisic_moments=False):
        """
        Read orbit library from file datfil/{fileroot}.dat.bz2'

        Parameters
        ----------
        fileroot : string
            this will probably be either 'orblib' or 'orblibbox'
        return_instrisic_moments: boolean
            whether to return_instrisic_moments of the orblib

        Returns
        -------
        if return_instrisic_moments is False, this returns a tuple of type
        (Histogram, array) where the orbit library LOSVDs are stored in the
        Histogram object, and the 3D density of the orbits are stored in the
        array object.
        if return_instrisic_moments is True, returns a tuple
        (array, list) where the array stores the intrinsic momenmts of the
        orblib and the list contains the bin edges of the 3D grid.

        """
        cur_dir = os.getcwd()
        os.chdir(self.mod_dir)
        # unzip orblib to a temporary file with ml value attached
        # ml value is needed to prevent different processes clashing
        ml = self.parset['ml']
        tmpfname = f'datfil/{fileroot}_{ml}.dat'
        subprocess.run(f'bunzip2 -c datfil/{fileroot}.dat.bz2 > {tmpfname}',
                       shell=True)
        # read the fortran file
        orblibf = FortranFile(tmpfname, 'r')
        # read size of orbit library
        # from integrator_setup_write, lines 506 - 5129:
        norb, _, _, _, ndith = orblibf.read_ints(np.int32)
        # from qgrid_setup_write, lines 2339-1350:
        quad_light_grid_sizes = orblibf.read_ints(np.int32)
        size_ql, size_qph, size_qth, size_qlr = quad_light_grid_sizes
        quad_lr = orblibf.read_reals(float)
        quad_lth = orblibf.read_reals(float)
        quad_lph = orblibf.read_reals(float)
        if return_instrisic_moments:
            instrisic_grid = [quad_lr, quad_lth, quad_lph]
        # from histogram_setup_write, lines 1917-1926:
        tmp = orblibf.read_record(np.int32, np.int32, float)
        nconstr = tmp[0][0] # = total number of apertures for ALL kinematics
                            #   and populations with own apertures
        # nvhist = tmp[1][0] # = (nvbins-1)/2 for histo of FIRST kinematic set
        # dvhist = tmp[2][0] # = delta_v in histogram for FIRST kinematic set
        # these nvhist and dvhist are for the first kinematic set only
        # however, orbit are stored N times where N = number of kinematic sets
        # histogram settings for other N-1 sets may be different from the first
        # these aren't stored in orblib.dat so must read from kinematics objects
        stars = self.system.get_unique_triaxial_visible_component()
        # For now, also calculate velocity histograms for population data...
        n_kins = len(stars.kinematic_data)
        pops = [p for p in stars.population_data if p.kin_aper is None]
        hist_widths = [k.hist_width for k in stars.kinematic_data]
        hist_widths += [p.hist_width for p in pops]
        hist_centers = [k.hist_center for k in stars.kinematic_data]
        hist_centers += [p.hist_center for p in pops]
        hist_bins = [k.hist_bins for k in stars.kinematic_data]
        hist_bins += [p.hist_bins for p in pops]
        self.logger.debug('Checking number of velocity bins...')
        error_msg = 'must have odd number of velocity bins for all ' \
                    'kinematics and populations'
        assert np.all(np.array(hist_bins) % 1==0), error_msg
        self.logger.debug('...checks ok.')
        n_apertures = [k.n_spatial_bins for k in stars.kinematic_data]
        n_apertures += [p.n_spatial_bins for p in pops]
        # get index linking  kinematic set to aperture
        # kin_idx_per_ap[i] = N <--> aperture i is from kinematic set N
        # NEW: Extended to kinematic and population data, this now reads:
        # get index linking kinematic and population sets with own apertures
        # to apertures
        # kinpop_idx_per_ap[i] = N <--> aperture i is from kinematic set N
        #   if N < n_kins or from population set N - n_kins if N >= n_kins
        kinpop_idx_per_ap = [np.zeros(n_apertures[i], dtype=int) + i
                             for i in range(n_kins + len(pops))]
        kinpop_idx_per_ap = np.concatenate(kinpop_idx_per_ap)
        kinpop_idx_per_ap = np.array(kinpop_idx_per_ap, dtype=int)
        # below we loop i_ap from 1-n_total_apertures but will need the index of
        # i_ap for the relevant kinematic set: we use `idx_ap_reset` to do this
        cum_n_apertures = np.cumsum(n_apertures)
        idx_ap_reset = np.concatenate(([0], cum_n_apertures[:-1]))
        # set up a list of arrays to hold the results
        tmp = zip(hist_bins,n_apertures)
        velhist0 = [np.zeros((norb, nv, na)) for (nv,na) in tmp]
        # Next read the histograms themselves.
        orbtypes = np.zeros((norb, ndith**3), dtype=int)
        density_3D = np.zeros((norb, size_qlr, size_qth, size_qph))
        if return_instrisic_moments:
            intrinsic_moms = np.zeros((norb, size_qlr, size_qth, size_qph, 16))
        for j in range(norb):
            _, _, _, _, _ = orblibf.read_ints(np.int32)
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
                kinpop_idx = kinpop_idx_per_ap[i_ap]
                i_ap0 = i_ap - idx_ap_reset[kinpop_idx]
                ivmin, ivmax = orblibf.read_ints(np.int32)
                if ivmin <= ivmax:
                    nv0 = (hist_bins[kinpop_idx]-1)/2
                    # ^--- this is an integer since hist_bins is odd
                    nv0 = int(nv0)
                    tmp = orblibf.read_reals(float)
                    velhist0[kinpop_idx][j, ivmin+nv0:ivmax+nv0+1, i_ap0] = tmp
            if return_instrisic_moments:
                intrinsic_moms[j] = quad_light
        orblibf.close()
        # remove temporary file
        os.remove(tmpfname)
        os.chdir(cur_dir)
        # For now, also calculate velocity histograms for population data...
        velhists = []
        for i in range(n_kins + len(pops)):
            center0 = hist_centers[i]
            width0 = hist_widths[i]
            bins0 = hist_bins[i]
            idx_center = (bins0-1)/2 # this is an integer since hist_bins is odd
            idx_center = int(idx_center)
            dvhist0 = width0/bins0
            vedg = np.arange(bins0+1) * dvhist0
            v = (vedg[1:] + vedg[:-1])/2.
            v_cent = v[idx_center]
            vedg -= v_cent
            vvv = dyn_kin.Histogram(xedg=vedg,
                                    y=velhist0[i],
                                    normalise=False)
            velhists += [vvv]
        if return_instrisic_moments:
            return intrinsic_moms, instrisic_grid
        else:
            return velhists, density_3D

    def duplicate_flip_and_interlace_orblib(self, orblib):
        """flip the tube orbits

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

    def duplicate_flip_and_interlace_intmoms(self, intmom):
        """equiv of `duplicate_flip_and_interlace_orblib` for intrinsic moments
        """
        new_shape = (intmom.shape[0]*2,) + intmom.shape[1:]
        new_intmom = np.zeros(new_shape)
        new_intmom[0::2] = intmom
        reversed_intmom = 1.* intmom # hack to make a copy
        # flip sign of...
        reversed_intmom[:,:,:,:,4] *= -1. # ... vx
        reversed_intmom[:,:,:,:,5] *= -1. # ... vy
        reversed_intmom[:,:,:,:,6] *= -1. # ... vz
        new_intmom[1::2, :] = reversed_intmom
        return new_intmom

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

    def read_losvd_histograms(self, kins=True, pops=False):
        """Read the orbit library

        Read box orbits and tube orbits, mirrors the latter, and combines.
        Rescales the velocity axis according to the ``ml`` value. Sets LOSVDs
        and 3D grid/aperture masses of the combined orbit library.

        Parameters
        ----------
        kins : Bool
            If True, return LOSVD histograms for the kinematics apertures.
        pops : Bool
            If True, return LOSVD histograms for the populations apertures.
            If both kins and pops are True, population data is returned
            following kinematic data. Note that if population and kinematic
            data are provided in one file (i.e. they share the same binning)
            then they have identical orbit libraries.

        Returns
        -------
        Sets the attributes:
            -   ``self.losvd_histograms``: a list, whose i'th entry is a
                ``dyn.kinematics.Histogram`` object holding the orbit lib LOSVDs
                binned for the i'th kinematic/population set
            -   ``self.intrinsic_masses``: 3D grid/intrinsic masses of orbit lib
            -   ``self.projected_masses``: aperture/proj. masses of orbit lib
            -   ``self.n_orbs``: number of orbits in the orbit library

        """
        if self.system.is_bar_disk_system():
            stars = self.system.get_unique_bar_component()
        else:
            stars = self.system.get_unique_triaxial_visible_component()
        n_kins = len(stars.kinematic_data)

        # TODO: check if this ordering is compatible with weights read in by
        # LegacyWeightSolver.read_weights
        tube_orblib, tube_density_3D = self.read_orbit_base('orblib')
        if not (kins or pops):
            txt = 'Specify kins or pops.'
            self.logger.error(txt)
            raise ValueError(txt)
        if kins and not pops:
            tube_orblib = tube_orblib[:n_kins]
        else:  # pops is True
            # collect all orblibs used by pops:
            kin_aper_in_pops = set(p.kin_aper
                                   for p in stars.population_data
                                   if p.kin_aper is not None)
            # kin orblibs + kin orblibs used by pops + pop-only orblibs:
            tube_orblib = \
                (tube_orblib[:n_kins] if kins else []) + \
                [o for i, o in enumerate(tube_orblib[:n_kins])
                   if i in kin_aper_in_pops] + \
                (tube_orblib[n_kins:] if len(tube_orblib) > n_kins else [])

        if not self.system.is_bar_disk_system():
            # tube orbits are mirrored/flipped and used twice
            tmp = []
            for tube_orblib0 in tube_orblib:
                tmp += [self.duplicate_flip_and_interlace_orblib(tube_orblib0)]
            tube_orblib = tmp
            tube_density_3D = np.repeat(tube_density_3D, 2, axis=0)

        # read box orbits
        box_orblib, box_density_3D = self.read_orbit_base('orblibbox')
        if kins and not pops:
            box_orblib = box_orblib[:n_kins]
        else:  # pops is True
            # kin orblibs + kin orblibs used by pops + pop-only orblibs:
            box_orblib = \
                (box_orblib[:n_kins] if kins else []) + \
                [o for i, o in enumerate(box_orblib[:n_kins])
                   if i in kin_aper_in_pops] + \
                (box_orblib[n_kins:] if len(box_orblib) > n_kins else [])

        # combine orblibs
        orblib = []
        for (t0, b0) in zip(tube_orblib, box_orblib):
            orblib0 = self.combine_orblibs(t0, b0)
            orblib += [orblib0]
        # combine density_3D arrays
        density_3D = np.vstack((tube_density_3D, box_density_3D))
        nkins = len(orblib)
        for i in range(nkins):
            orblib[i].scale_x_values(self.velocity_scaling_factor)
        self.losvd_histograms = orblib
        self.intrinsic_masses = density_3D
        self.n_orbs = self.losvd_histograms[0].y.shape[0] if nkins > 0 else 0
        proj_mass = [np.sum(self.losvd_histograms[i].y,1) for i in range(nkins)]
        self.projected_masses = proj_mass

    def read_orbit_intrinsic_moments(self, cache=True):
        """Read the intrinsic moments of the orbit library.

        Moments stored in 3D grid over spherical co-ords (r,theta,phi). This
        function reads the data from files, formats them correctly, and coverts
        to physical units.

        Parameters
        ----------
        cache : bool, optional
            If True, the intrinsic moments and the bin edges are cached
            between calls to this method. The cache files are stored in the
            orbit library's datfil/ directory. The default is True.

        Returns
        -------
        (array, list)
            array shape = (n_orb, nr, nth, nph, 16). Final dimension indexes
            over: density,x,y,z,vx,vy,vz,vx^2,vy^2,vz^2,vx*vy,vy*vz,vz*vx, and
            the final three indices (13,14,15) are some type of orbit
            classification (not understood - recommend not to use!). The list
            contains grid bin edges over spherical (r, theta, phi).

        """
        mom_file = self.mod_dir + 'datfil/intmoms.npz'
        grid_file = self.mod_dir + 'datfil/int_grid.npz'
        mom_file_exists = os.path.isfile(mom_file)
        grid_file_exists = os.path.isfile(grid_file)
        if cache and mom_file_exists and grid_file_exists:
            with np.load(mom_file, allow_pickle=False) as data:
                intmoms = data['intmoms']
            with np.load(grid_file, allow_pickle=False) as data:
                int_grid = [data[x] for x in ('r', 'th', 'ph')]
        else:
            intmom_tubes, int_grid = self.read_orbit_base(
                'orblib',
                return_instrisic_moments=True)
            intmom_tubes = \
                self.duplicate_flip_and_interlace_intmoms(intmom_tubes)
            intmom_boxes, _ = self.read_orbit_base(
                'orblibbox',
                return_instrisic_moments=True)
            intmoms = np.concatenate((intmom_tubes, intmom_boxes), 0)
            velscale = self.velocity_scaling_factor
            conversion_factor = self.system.distMPc * 1e6 * \
                                np.tan(np.pi/648000.0) * PARSEC_KM
            intmoms[:,:,:,:,1] /= conversion_factor # kpc -> arcsec
            intmoms[:,:,:,:,2] /= conversion_factor
            intmoms[:,:,:,:,3] /= conversion_factor
            intmoms[:,:,:,:,4] *= velscale # for velocity stretching due to M/L
            intmoms[:,:,:,:,5] *= velscale
            intmoms[:,:,:,:,6] *= velscale
            intmoms[:,:,:,:,7] *= velscale**2.
            intmoms[:,:,:,:,8] *= velscale**2.
            intmoms[:,:,:,:,9] *= velscale**2.
            intmoms[:,:,:,:,10] *= velscale**2.
            intmoms[:,:,:,:,11] *= velscale**2.
            intmoms[:,:,:,:,12] *= velscale**2.
            int_grid[0] /= conversion_factor # kpc -> arcsec
            if cache:
                if mom_file_exists:
                    os.remove(mom_file)
                if grid_file_exists:
                    os.remove(grid_file)
                np.savez_compressed(mom_file, intmoms=intmoms)
                np.savez(grid_file,
                         r=int_grid[0],
                         th=int_grid[1],
                         ph=int_grid[2])
        return intmoms, int_grid

    def read_orbit_property_file_base(self, file, ncol, nrow):
        """Base method to read in ``*orbclass.out`` files

        ...which hold the information of all the orbits stored in the orbit
        library. The number of orbits is ``nE * nI2 * nI3 * dithering^3``,
        ``ncol = dithering^3`` and ``nrow = nE * nI2 * nI3``.
        For each orbit, the time averaged values are stored:
        ``lx, ly ,lz, r = sum(sqrt( average(r^2) ))``,
        ``Vrms^2 = average(vx^2 + vy^2 + vz^2 + 2vx*vy + 2vxvz + 2vxvy)``.
        The files were originally stored by the fortran code ``orblib_f.f90``,
        ``integrator_find_orbtype``.

        *** Do not try to replace this with numpy (e.g.** ``np.genfromtext`` **)
        since the output files are sometimes written in non-standard arrays
        (in a seemingly system dependent way). ***

        """
        data=[]
        lines = [line.rstrip('\n').split() for line in open(file)]
        i = 0
        while i < len(lines):
            for x in lines[i]:
                data.append(np.double(x))
            i += 1
        data=np.array(data)
        if len(data) != 5*ncol*nrow:
            txt = f'{file} length mismatch - found {len(data)} entries, ' \
                  f'expected: {5*ncol*nrow}. Correct configuration file used?'
            self.logger.error(txt)
            raise ValueError(txt)
        data=data.reshape((5,ncol,nrow), order='F')
        return data

    def read_orbit_property_file(self):
        """Read the ``*orbclass.out`` files

        These files contain time-averaged properties of individual orbits within
        a bundle. Results are stored in ``self.orb_properties``, an astropy
        table with columns ``[r, Vrms, L, Lx, Ly, Lz, lmd, lmd_x, lmd_y,
        lmd_z]`` i.e.
        the time-averaged value of (i) radius ``r``, (ii) ``V_rms =
        (vx^2 + vy^2 + vz^2 + 2(vxvy + vxvz + vyvz))``, (iii) total angular
        momentum ``L``, (iv) ``Lx``, (v) ``Ly``, (vi) ``Lz``, and
        circularities (vii) ``lmd = L/V_rms``, (ix) ``lmd_x=Lx/V_rms``,
        (x) ``lmd_y``, (xi) ``lmd_z``. Each column
        of the table is an array of shape ``[N_bundle, N_orbit]``.

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

    def classification_diagnostic_plot(self, dl=1.0e17):
        """Find threshold angular momentum ``dL`` for use in orbit classification

        Parameters
        ----------
        dl : float, threshold angular momentum, default = 1e17

        Returns
        -------
        dl : float
            The threshold angular momentum ``dL``.
        """
        orb_properties = self.orb_properties
        kw_hist = {'range':(-1.0e18, 1.0e18), 'bins':51, 'alpha':0.3}
        _ = plt.hist(np.ravel(orb_properties['Lx'].value), label='$L_x$', **kw_hist)
        _ = plt.hist(np.ravel(orb_properties['Ly'].value), label='$L_y$', **kw_hist)
        _ = plt.hist(np.ravel(orb_properties['Lz'].value), label='$L_z$', **kw_hist)
        plt.axvline(dl, ls=':', color='k')
        plt.axvline(-dl, ls=':', color='k')
        plt.gca().set_xlabel('Angular momenta $L_{x/y/z}$ [km/s/kpc]')
        plt.gca().set_ylabel('pdf')
        plt.gca().set_title('Distribtuion of angular momenta of orbits in orblib')
        plt.gca().legend()
        return dl

    def classify_orbits(self, make_diagnostic_plots=False, dL=1.0e17, force_lambda_z=False):
        """ Perform orbit classification.

        Orbits are classified in 3D angular momentum space ``(Lx, Ly, Lz)``
        according to a threshold angular momdentum ``dL``. Ideally, we would
        classify as follows:

        - ``|Lx|,|Ly|,|Lz|<dL`` --> box
        - ``|Lx|>dL`` & ``|Ly|,|Lz|<dL`` --> x-axis tube
        - ``|Ly|>dL`` & ``|Lx|,|Lz|<dL`` --> y-axis tube (theoretically unstable)
        - ``|Lz|>dL`` & ``|Lx|,|Ly|<dL`` --> z-axis tube

        but this "exact" classification leaves many orbits unclassified.
        Instead we classify tube orbits using the less strict criteria:

        - ``|Lx|,|Ly|,|Lz|<dL`` --> box
        - not a box & ``|Lx|>|Ly|`` & ``|Lx|>|Lz|`` --> x-axis tube "-ish"
        - not a box & ``|Ly|>|Lz|`` & ``|Ly|>|Lz|`` --> y-axis tube "-ish"
        - not a box & ``|Lz|>|Lx|`` & ``|Lz|>|Ly|`` --> z-axis tube "-ish"

        The method logs the fraction of all orbits in each classification, and
        the fraction of each type which are "exact" according to the stricter
        criteria. The results saved in ``self.orb_classification`` are for the
        less strict criteria.

        Parameters
        ----------
        make_diagnostic_plots : bool, optional
            whether to make diagnostic plot in
            ``find_threshold_angular_momentum``, by default ``False``.

        dL : float, threshold angular momentum, default = 1e17

        Returns
        -------
        None
            sets result to ``self.orb_classification``, a dictionary containing
            5 keys ``bool_{box,xtish,ytish,ztish,other}`` where e.g.
            ``bool_box`` is a boolean array with shape ``[N_bundles, N_orbits]``
            set to ``True`` at the location of box orbits.
            ``bool_other`` should be ``False`` everywhere.

        """
        orb_properties = self.orb_properties
        dl = 1.*dL
        if make_diagnostic_plots:
            self.classification_diagnostic_plot(dl=dl)
        # find box orbits
        bool_box = (
            (np.abs(orb_properties['Lx'].value) <= dl) &
            (np.abs(orb_properties['Ly'].value) <= dl) &
            (np.abs(orb_properties['Lz'].value) <= dl)
        )
        idx_box = np.where(bool_box)
        # find "true" tube orbits i.e. with exactly one component of L =/= 0
        bool_xtube = (
            (np.abs(orb_properties['Lx'].value) > dl) &
            (np.abs(orb_properties['Ly'].value) <= dl) &
            (np.abs(orb_properties['Lz'].value) <= dl)
        )
        bool_ytube = (
            (np.abs(orb_properties['Lx'].value) <= dl) &
            (np.abs(orb_properties['Ly'].value) > dl) &
            (np.abs(orb_properties['Lz'].value) <= dl)
        )
        bool_ztube = (
            (np.abs(orb_properties['Lx'].value) <= dl) &
            (np.abs(orb_properties['Ly'].value) <= dl) &
            (np.abs(orb_properties['Lz'].value) > dl)
        )
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
        if force_lambda_z:
            bool_xtish = np.full_like(bool_xtish, False)
            bool_ytish = np.full_like(bool_ytish, False)
            bool_box = np.full_like(bool_box, False)
            bool_ztish = np.full_like(bool_ztish, True)
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
        frac = n_xt_exact/n_xtish if n_xtish > 0 else 0
        self.logger.info(f'    - {percent(frac)} of x-tubes')
        frac = n_yt_exact/n_ytish if n_ytish > 0 else 0
        self.logger.info(f'    - {percent(frac)} of y-tubes')
        frac = n_zt_exact/n_ztish if n_ztish > 0 else 0
        self.logger.info(f'    - {percent(frac)} of z-tubes')
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

    def get_projection_tensor(self, minr=None, maxr=None, nr=50, nl=61, force_lambda_z=False, dL=1e17):
        projection_tensor_pars = {'minr':minr,
                                  'maxr':maxr,
                                  'nr':nr,
                                  'nl':nl,
                                  'dL':dL,
                                  'force_lambda_z':force_lambda_z}
        self.projection_tensor_pars = projection_tensor_pars
        # otherwise, continue...
        if hasattr(self, 'orb_properties') == False:
            self.read_orbit_property_file()
        orb_properties = self.orb_properties
        self.classify_orbits(dL=dL, force_lambda_z=force_lambda_z)
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
            projection00 = np.mean(projection00, axis=1)
            projection += [projection00]
        projection = np.stack(projection)
        projection = np.moveaxis(projection, 1, 3)
        self.projection_tensor = projection

    def get_model_intrinsic_moment_constructor(self):
        """Get a function to constrcut the model's intrinsic moments in 3D grid

        Returns a function and a list. The function takes weights and returns
        intrinsic model moments in a 3D grid. The list contains the bin edges
        of the 3D grid. Example usage:

        ```
        moment_constructor, bin_edges = orblib.get_model_intrinsic_moment_constructor()
        moments = moment_constructor(weights)
        ```

        Returns
        -------
        (callable, list)
            the callable = function which takes weights and returns 3D moments.
            The list contains grid bin edges over spherical (r, theta, phi).
            Moments returned by the callable are stored in a grid of size
            (nr, nth, nph, 13). Final dimension indexes over: density,x,y,z,vx,
            vy,vz,vx^2,vy^2,vz^2,vx*vy,vy*vz,vz*vx. Density is normalised to 1,
            spatial moments in arcseconds, velocities in km/s.

        """
        intmoms, int_grid = self.read_orbit_intrinsic_moments()
        density = intmoms[:,:,:,:,0]
        kinmoms = intmoms[:,:,:,:,1:13]
        def model_intrinsic_moment_constructor(weights):
            mod_density = (density.T * weights).T
            # normalise density so model sums to 1 in each (r,th,ph) bin
            mod_dens_nrm = mod_density/np.sum(mod_density, 0)
            mod_kinmoms = np.einsum('ijklm,ijkl->jklm', kinmoms, mod_dens_nrm)
            # normalise density that model sums over (r,th,ph) bins to give 1
            mod_density = np.sum(mod_density, 0)
            mod_density /= np.sum(mod_density)
            mod_moments = np.concatenate(
                (mod_density[...,np.newaxis], mod_kinmoms),
                axis=-1)
            return mod_moments
        return model_intrinsic_moment_constructor, int_grid

# end
