import os
import pathlib
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
from dynamite import orblib_chunks
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
        self.config = config
        self.mod_dir = mod_dir
        self.parset = parset
        self.system = config.system
        self.settings = config.settings.orblib_settings
        self.legacy_directory = config.settings.legacy_settings['directory']
        self.in_dir = config.settings.io_settings['input_directory']
        weight_solver = config.settings.weight_solver_settings['type']
        self.LegacyWeightSolver = weight_solver == 'LegacyWeightSolver'
        self.orblibs_in_parallel = \
            config.settings.multiprocessing_settings['orblibs_in_parallel']
        self.n_chunks = \
            config.settings.multiprocessing_settings['orblib_chunks']
        self.stars = self.system.get_unique_triaxial_visible_component()
        self.n_hist1d = len([k for k in self.stars.kinematic_data
                             if not isinstance(k, dyn_kin.ProperMotions)])
        self.n_hist2d = len([k for k in self.stars.kinematic_data
                             if isinstance(k, dyn_kin.ProperMotions)])
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
            - orblib_qgrid.dat.bz2          (zipped tube orbit library)
            - orblib_losvd_hist.dat.bz2(*)  (zipped losvd 1d histograms)
            - orblib_pops.dat.bz2(*)        (zipped populations 0d hists)
            - orblib_pm_hist.dat.bz2(*)     (zipped proper motions 2d histograms)
            - orblib.dat_orbclass.out       (orbit classification for tube orbs)
            - orblibbox_qgrid.dat.bz2(*)    (zipped box orbit library)
            - orblibbox_losvd_hist.dat.bz2(*) (zipped losvd 1d histograms)
            - orblibbox_pops.dat.bz2(*)     (zipped populations 0d hists)
            - orblibbox_pm_hist.dat.bz2(*)  (zipped proper motions 2d histograms)
            - orblibbox.dat_orbclass.out    (orbit classification for box orbs)
            - LegacyWeightSolver only:
                - mass_aper.dat             (MGE masses in apertures)
                - mass_qgrid.dat            (MGE masses in 3D grid)
                - mass_radmass.dat          (MGE masses in radial bins)
            - Not for LegacyWeightSolver:
                - mass_qgrid.ecsv           (MGE masses in 3D grid)
                - mass_radmass.ecsv         (MGE masses in radial bins)
            - + up to 8 log and status files
        (*) only created if there are kinematic datasets with losvd 1d
            histograms, populations datasets, or proper motion datasets, resp.

        """
        # check whether orbit library was calculated already
        if not os.path.isfile(self.mod_dir + 'datfil/tube_box_done'):
            # prepare the fortran input files for orblib
            self.create_fortran_input_orblib(self.mod_dir+'infil/')
            stars = self.stars
            # create the kinematics and populations input files for each
            # kinematic dataset and population dataset with own apertures
            pops = [p for p in stars.population_data if p.kin_aper is None]
            for data_set in stars.kinematic_data + pops:
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
            if self.n_chunks > 1 and self.can_chunk_orbits():
                self.get_orbit_library_chunked(self.n_chunks)
            elif self.orblibs_in_parallel:
                self.get_orbit_library_par()
            else:
                self.get_orbit_library()
            if not self.LegacyWeightSolver:
                # Calculate the orblib's parset's observed intrinsic masses
                model=self.config.all_models.get_model_from_parset(self.parset)
                if self.system.is_bar_disk_system():
                    mge = stars.mge_lum_tot
                    len_mge_bulge = len(stars.mge_lum.data)
                    # intrinsic masses
                    _ = mge.get_intrinsic_masses(model,
                                                 len_mge_bulge=len_mge_bulge,
                                                 parallel=False)
                else:
                    mge = stars.mge_lum
                    # intrinsic masses
                    _ = mge.get_intrinsic_masses(model, parallel=False)
            tube_done = os.path.isfile(self.mod_dir + 'datfil/tube_done')
            box_done = os.path.isfile(self.mod_dir + 'datfil/box_done')
            if tube_done and box_done:
                pathlib.Path(self.mod_dir + 'datfil/tube_box_done').touch()


    def create_fortran_input_orblib(self, path):
        """write input files for Fortran orbit library programs

        Parameters
        ----------
        path : string
            path of the model's ``infil`` directory

        Returns
        -------
        Creates the following files in the ``infil`` directory:
            - parameters_pot.in
            - parameters_lum.in
            - orblib.in
            - orblibbox.in
            - triaxmass.in (LegacyWeightSolver only)
            - triaxmassbin.in (LegacyWeightSolver only)

        """
        #---------------------------------------------
        #write parameters_pot.in and parameters_lum.in
        #---------------------------------------------
        stars = self.stars
        bh = self.system.get_component_from_class(physys.Plummer)
        # used to derive the viewing angles
        if self.system.is_bar_disk_system():
            theta = self.parset[f'theta-{stars.name}']
            phi = self.parset[f'phi-{stars.name}']
            psi = self.parset[f'psi-{stars.name}']
        else:
            q = self.parset[f'q-{stars.name}']
            p = self.parset[f'p-{stars.name}']
            u = self.parset[f'u-{stars.name}']
            theta, psi, phi = stars.triax_pqu2tpp(p,q,u)
        # get dark halo
        dh = self.system.get_all_dark_non_plummer_components()
        if len(dh) > 1:
            txt = 'Zero or one non-plummer dark component should be ' \
                  f' present, not {len(dh)}.'
            self.logger.error(txt)
            raise ValueError(txt)
        if len(dh) > 0:
            dh = dh[0]  # extract the one and only dm component

            if isinstance(dh, physys.NFW_m200_c):
                #fix c via m200_c relation, for legacy Fortran it is still NFW
                dm_specs, dm_par_vals = dh.get_dh_legacy_strings(self.parset,
                                                                self.system)
            else:
                dm_specs, dm_par_vals = dh.get_dh_legacy_strings(self.parset)
        else:
            dm_specs = '0 0'
            dm_par_vals = ''

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
        text += f"{settngs['quad_nr']}\n"
        text += f"{settngs['quad_nth']}\n"
        text += f"{settngs['quad_nph']}\n"
        text += f"{dm_specs}\n"
        text += f"{dm_par_vals}\n"
        if self.system.is_bar_disk_system():
            len_disk_pot = len(stars.disk_pot.data)
            header_string_pot = str(len_mge_pot + len_disk_pot) + " 1 " + str(len_mge_pot) + " " + str(len_disk_pot)
            len_disk_lum = len(stars.disk_lum.data)
            header_string_lum = str(len_mge_lum + len_disk_lum) + " 1 " + str(len_mge_lum) + " " + str(len_disk_lum)
            text += f"{self.parset['omega']}\n"
            mge_pot = stars.mge_pot + stars.disk_pot
            mge_lum = stars.mge_lum + stars.disk_lum
        else:
            header_string_pot = str(len_mge_pot)
            header_string_lum = str(len_mge_lum)
            mge_pot = stars.mge_pot
            mge_lum = stars.mge_lum
        text += f"{self.system.H*1e-6}"

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
        n_psf_kin = len(stars.kinematic_data)  # write data for all kinematics
        psf_pop_idx = [i for i, p in enumerate(stars.population_data)
                       if p.kin_aper is None]  # pops with their own apertures
        def write_orblib_dot_in(box=False, start=None, number=None, tag=''):
            """Write one orbit library input file.

            Parameters
            ----------
            box : bool, optional
                write the box orbit file rather than the tube one.
            start, number : int, optional
                first orbit and orbit count, overriding the ``starting_orbit``
                and ``number_orbits`` settings. Used to give each chunk its
                own slice of the library. The defaults are None, meaning use
                the settings.
            tag : string, optional
                suffix for this file and for the orbit library it writes, so
                that chunks integrated concurrently do not collide. The
                default is '', i.e. the usual unchunked file names.

            """
            tab = '\t\t\t\t\t\t\t\t'
            fileroot = 'orblibbox' if box else 'orblib'
            f = open(f'{path}{fileroot}{tag}.in', 'w')
            label = '[random seed for orbit integration]'
            f.write(f"{self.settings['random_seed']}{tab}{label}\n")
            f.write('infil/parameters_pot.in\n')
            if box:
                f.write('datfil/beginbox.dat\n')
            else:
                f.write('datfil/begin.dat\n')
            label = '[number of orbital periods to integrate]'
            line = f"{self.settings['orbital_periods']}{tab}{label}\n"
            f.write(line)
            label = '[points to sample for each orbit in the merid. plane]'
            line = f"{self.settings['sampling']}{tab}{label}\n"
            f.write(line)
            label = '[starting orbit]'
            s = self.settings['starting_orbit'] if start is None else start
            line = f"{s}{tab}{label}\n"
            f.write(line)
            label = '[orbits to integrate; -1 --> all orbits]'
            n = self.settings['number_orbits'] if number is None else number
            line = f"{n}{tab}{label}\n"
            f.write(line)
            label = '[accuracy]'
            line = f"{self.settings['accuracy']}{tab}{label}\n"
            f.write(line)
            label = '[number of psfs of the kinematic data + proper motions + population data]'
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
            label = '[# of apertures of the kinematic data + proper motion data + population data]'
            line = f'{n_aperture_kin + n_aperture_pop}{tab}{label}\n'
            f.write(line)
            for i in np.arange(n_aperture_kin):
                kin_i = stars.kinematic_data[i]
                # aperturefile cant have a label since fortran reads whole line
                f.write(f'"infil/{kin_i.aperturefile}"\n')
                label = f'[use kinematics psf {i+1} for {kin_i.name}]'
                line = f'{i+1}{tab}{label}\n'
                f.write(line)
                label = f'[histogram dimension for apert {1+i} ({kin_i.name})]'
                ap_hist_dim = 2 if isinstance(kin_i, dyn_kin.ProperMotions) \
                                else 1
                f.write(f'{ap_hist_dim}{tab}{label}\n')
            for i0, i in enumerate(psf_pop_idx):  # i0 counts aperture sets
                pop_i = stars.population_data[i]
                # aperturefile cant have a label since fortran reads whole line
                f.write(f'"infil/{pop_i.aperturefile}"\n')
                label = f'[use populations psf {i+1} = ' \
                        f'total psf {i0+n_psf_kin+1} for {pop_i.name}]'
                line = f'{i0+n_psf_kin+1}{tab}{label}\n'
                f.write(line)
                label = f'[histogram dimension for psf {i0+n_psf_kin+1} ' \
                        f'({pop_i.name})]'
                ap_hist_dim = 0  # 0d histograms (scalars) for populations data
                f.write(f'{ap_hist_dim}{tab}{label}\n')
            for kin_i in stars.kinematic_data:
                txt = '[vhist width, center and nbins for kinematics ' \
                      f'{kin_i.name}]'
                if isinstance(kin_i, dyn_kin.ProperMotions):
                    for i in range(2):
                        label = f'dim {i + 1} of 2 {txt}'
                        w = kin_i.hist_width[i]
                        c = kin_i.hist_center[i]
                        b = kin_i.hist_bins[i]
                        line = f'{w} {c} {b} {tab[:-2]}{label}\n'
                        f.write(line)
                else:
                    label = f'dim 1 of 1 {txt}'
                    w, c = kin_i.hist_width, kin_i.hist_center
                    b = kin_i.hist_bins
                    line = f'{w} {c} {b}{tab[:-2]}{label}\n'
                    f.write(line)
            for i in psf_pop_idx:
                pop_i = stars.population_data[i]
                label = '[vhist width, center, nbins are trivial for pops ' \
                        f'{pop_i.name}]'
                w, c, b = 1., 0., 1
                line = f'{w} {c} {b}{tab[:-1]}{label}\n'
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
            o_file = f'datfil/{fileroot}{tag}'
            f_name = f'"{o_file}_qgrid.dat"'
            f.write(f'{f_name}{tab[:-4] if len(f_name) >= 32 else tab[:-3]}'
                    '[orbit qgrid file]\n')
            if len(psf_pop_idx) > 0:
                f_name = f'"{o_file}_pops.dat"'
                f.write(f'{f_name}{tab[:-4] if len(f_name) >= 32 else tab[:-3]}'
                        '[pops \'0d hist\' file]\n')
            if self.n_hist1d > 0:
                f_name = f'"{o_file}_losvd_hist.dat"'
                f.write(f'{f_name}{tab[:-4] if len(f_name) >= 32 else tab[:-3]}'
                        '[orbit losvd 1d hist file]\n')
            if self.n_hist2d > 0:
                f_name = f'"{o_file}_pm_hist.dat"'
                f.write(f'{f_name}{tab[:-4] if len(f_name) >= 32 else tab[:-3]}'
                        '[orbit proper motions 2d hist file]\n')
            f_name = f'"{o_file}.dat_orbclass.out"'
            f.write(f'{f_name}{tab[:-4] if len(f_name) >= 32 else tab[:-3]}'
                    '[orbit classification file]\n')
            f.close()
        write_orblib_dot_in(box=False)
        write_orblib_dot_in(box=True)
        # when integrating in chunks, each process needs its own input file
        # naming its own orbit range and its own output files
        if self.n_chunks > 1 and self.can_chunk_orbits():
            for start, number, tag in self.orbit_chunk_bounds(self.n_chunks):
                write_orblib_dot_in(box=False, start=start, number=number,
                                    tag=tag)
                write_orblib_dot_in(box=True, start=start, number=number,
                                    tag=tag)
        if self.LegacyWeightSolver:
            #--------------------------------------------
            #write triaxmass.in (LegacyWeightSolver only)
            #--------------------------------------------
            text = 'infil/parameters_lum.in\n' + \
                   'datfil/orblib_qgrid.dat\n' + \
                   'datfil/mass_radmass.dat\n' + \
                   'datfil/mass_qgrid.dat'
            triaxmass_file = open(path+'triaxmass.in',"w")
            triaxmass_file.write(text)
            triaxmass_file.close()
            #-----------------------------------------------
            #write triaxmassbin.in (LegacyWeightSolver only)
            #-----------------------------------------------
            tab = '\t\t\t\t\t\t\t\t'
            f = open(path + 'triaxmassbin.in', 'w')
            f.write('infil/parameters_lum.in\n')
            f.write(f'{n_psf_kin}{tab}[# of kinematics apertures]\n')
            for i in range(n_psf_kin):  # note: no pops here
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
            f.write('"datfil/mass_aper.dat"')
            f.close()

    def get_orbit_ics(self):
        """Run the Fortran executable to calculate orbit ICs
        """
        cur_dir = os.getcwd()
        os.chdir(self.mod_dir)
        self.logger.info(f'Calculating initial conditions for {self.mod_dir}.')
        cmd = self.legacy_directory
        bar = '_bar' if self.system.is_bar_disk_system() else ''
        cmd += f'/orbitstart{bar} < infil/orbstart.in >> datfil/orbstart.log'
        p = subprocess.run(cmd,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,
                           shell=True)
        os.chdir(cur_dir)
        log_file = f'Logfile: {self.mod_dir}datfil/orbstart.log.'
        if not p.stdout.decode("UTF-8"):
            self.logger.info(f'...done - orbitstart{bar} exit code '
                             f'{p.returncode}. {log_file}')
        else:
            text = f'...failed! orbitstart{bar} exit code {p.returncode}. ' \
                   f'Message: {p.stdout.decode("UTF-8")}. {log_file}'
            if p.returncode == 127: # command not found
                text += 'Check DYNAMITE legacy_fortran executables.'
                self.logger.error(text)
                raise FileNotFoundError(text)
            else:
                text += f'{log_file} Be wary: DYNAMITE may crash...'
                self.logger.warning(text)
                raise RuntimeError(text)

    def get_orbit_library_chunked(self, n_chunks):
        """Integrate the orbit library in chunks, then merge the results.

        Splits each orbit family's starting points into ``n_chunks``
        consecutive ranges, integrates all of them concurrently, merges each
        family's output back into the usual single files and compresses them.

        The merged library is **bit-identical** to the one produced by the
        unchunked path: the orbit range only selects which starting points a
        process integrates, and each orbit's result depends on its own index
        rather than on how many orbits preceded it.

        Parameters
        ----------
        n_chunks : int
            number of chunks per orbit family.

        Raises
        ------
        FileNotFoundError
            if a Fortran executable or an expected chunk file is missing.
        RuntimeError
            if the integration script reports an error.

        """
        cur_dir = os.getcwd()
        os.chdir(self.mod_dir)
        try:
            bounds = self.orbit_chunk_bounds(n_chunks)
            cmdstr = self.write_executable_for_integrate_orbits_chunked(bounds)
            self.logger.info(
                f'Integrating orbit library for {self.mod_dir} in '
                f'{len(bounds)} chunk(s) per family '
                f'({2 * len(bounds)} processes).')
            p = subprocess.run('bash ' + cmdstr,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT,
                               shell=True)
            if p.stdout.decode('UTF-8'):
                text = f'...failed! {cmdstr} exit code {p.returncode}. ' \
                       f'Message: {p.stdout.decode("UTF-8")}'
                if p.returncode == 127:
                    text += 'Check DYNAMITE legacy_fortran executables.'
                    self.logger.error(text)
                    raise FileNotFoundError(text)
                self.logger.warning(text + ' Be wary: DYNAMITE may crash...')
                raise RuntimeError(text)
            # merge before compressing: the chunks are only meaningful together
            tags = [tag for _, _, tag in bounds]
            n_orbits = self.n_orbit_starting_points() \
                * self.settings['dithering'] ** 3
            for fileroot in ('orblib', 'orblibbox'):
                # Models sharing an orbit library may be evaluated
                # concurrently, in which case more than one of them integrates
                # into this directory. If another has already merged and
                # compressed this family, its chunks are gone and there is
                # nothing left to do.
                done = [f'datfil/{fileroot}_{kind}.dat.bz2'
                        for kind in orblib_chunks.HEADER_RECORDS]
                if all(os.path.isfile(f) for f in done):
                    self.logger.debug(
                        f'{self.mod_dir}: {fileroot} already merged by a '
                        'concurrent model, discarding this process\' chunks.')
                    orblib_chunks.remove_chunks('datfil', fileroot, tags)
                    continue
                # Every intermediate is per-process; only the final .bz2 has a
                # shared name, and it appears atomically via mv. Concurrent
                # models therefore overwrite each other with identical content,
                # exactly as the unchunked path does, instead of colliding
                # halfway through.
                out_tag = f'.p{os.getpid()}'
                merged = orblib_chunks.merge_chunks(
                    'datfil', fileroot, tags, n_orbits, out_tag=out_tag)
                for f_name in merged:
                    final = f_name[:-len(out_tag)] + '.bz2'
                    subprocess.run(
                        f'bzip2 -c {f_name} > {f_name}.bz2 '
                        f'&& mv {f_name}.bz2 {final} && rm -f {f_name}',
                        shell=True, check=True)
                self.logger.debug(f'{self.mod_dir}: merged {len(tags)} '
                                  f'{fileroot} chunks.')
            self.logger.info(f'...done - {cmdstr} exit code {p.returncode}. '
                             f'Logfiles: {self.mod_dir}datfil/orblib.log, '
                             f'{self.mod_dir}datfil/orblibbox.log.')
        finally:
            os.chdir(cur_dir)

    def get_orbit_library_par(self):
        """Execute the bash script to calculate orbit libraries in parallel
        """
        # move to model directory
        cur_dir = os.getcwd()
        os.chdir(self.mod_dir)
        cmdstr = self.write_executable_for_integrate_orbits_par()
        self.logger.info('Integrating orbit library tube and box orbits '
                         f'for {self.mod_dir}.')
        # p = subprocess.call('bash '+cmdstr_tube, shell=True)
        p = subprocess.run('bash '+cmdstr,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,
                           shell=True)
        # move back to original directory
        os.chdir(cur_dir)
        log_files = f'Logfiles: {self.mod_dir}datfil/orblib.log, ' \
                    f'{self.mod_dir}datfil/orblibbox.log'
        if self.LegacyWeightSolver:
            log_files += f', {self.mod_dir}datfil/triaxmassbin.log'
            log_files += f', {self.mod_dir}datfil/triaxmass.log.'
        else:
            log_files += '.'
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

    def get_orbit_library(self):
        """Execute the bash script to calculate orbit libraries: first tube,
        then box orbits
        """
        # move to model directory
        cur_dir = os.getcwd()
        os.chdir(self.mod_dir)
        cmdstr_tube, cmdstr_box = self.write_executable_for_integrate_orbits()
        self.logger.info('Integrating orbit library tube orbits '
                         f'for {self.mod_dir}.')
        # p = subprocess.call('bash '+cmdstr_tube, shell=True)
        p = subprocess.run('bash '+cmdstr_tube,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,
                           shell=True)
        # move back to original directory
        os.chdir(cur_dir)
        log_files = f'Logfiles: {self.mod_dir}datfil/orblib.log'
        if self.LegacyWeightSolver:
            log_files += f', {self.mod_dir}datfil/triaxmassbin.log'
            log_files += f', {self.mod_dir}datfil/triaxmass.log.'
        else:
            log_files += '.'
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
        os.chdir(self.mod_dir)
        self.logger.info('Integrating orbit library box orbits '
                         f'for {self.mod_dir}.')
        # p = subprocess.call('bash '+cmdstr_box, shell=True)
        p = subprocess.run('bash '+cmdstr_box,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,
                           shell=True)
        # move back to original directory
        os.chdir(cur_dir)
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

    def n_orbit_starting_points(self):
        """Number of orbit starting points in the library.

        This is what the Fortran's ``starting orbit`` and ``orbits to
        integrate`` input lines count, and therefore what a chunk range refers
        to. Each starting point expands into ``dithering^3`` integrated orbits.

        Returns
        -------
        int
            ``nE * nI2 * nI3 / dithering^3``

        """
        s = self.settings
        return (s['nE'] * s['nI2'] * s['nI3']) // s['dithering'] ** 3

    def can_chunk_orbits(self):
        """Whether this orbit library may be integrated in chunks.

        Chunking merges the ``qgrid`` and ``losvd_hist`` output streams. The
        proper motion (2d histogram) and populations streams are not merged, so
        libraries producing them are integrated in one piece regardless of the
        ``orblib_chunks`` setting.

        Returns
        -------
        bool
            True if chunked integration is supported for this library.

        """
        has_pops = any(p.kin_aper is None for p in self.stars.population_data)
        return self.n_hist2d == 0 and not has_pops

    def orbit_chunk_bounds(self, n_chunks):
        """Divide the orbit starting points into consecutive chunks.

        The remainder is spread one per chunk rather than landing entirely on
        the last one, so chunk sizes differ by at most one orbit.

        Parameters
        ----------
        n_chunks : int
            number of chunks. Clamped to the number of starting points.

        Returns
        -------
        list of tuple
            one ``(start, number, tag)`` per chunk, in ascending orbit order.
            ``start`` is 1-based, matching the Fortran input file, and ``tag``
            is the suffix distinguishing this chunk's files.

        Notes
        -----
        The tag carries the process id. Models sharing an orbit library may be
        evaluated concurrently, and each such process integrates into the same
        ``datfil`` directory; without the process id they would write over each
        other's chunks, and one could merge another's half-written files. The
        merged output does not carry the tag, so concurrent processes still
        produce the same final library, exactly as the unchunked path does.

        """
        n_orbits = self.n_orbit_starting_points()
        n_chunks = max(1, min(int(n_chunks), n_orbits))
        base, rem = divmod(n_orbits, n_chunks)
        bounds, start = [], 1
        for k in range(n_chunks):
            number = base + (1 if k < rem else 0)
            bounds.append((start, number, f'_p{os.getpid()}_c{k}'))
            start += number
        return bounds

    def write_executable_for_integrate_orbits_chunked(self, bounds):
        """Write the bash script integrating both orbit families in chunks.

        Every chunk of both families runs concurrently, giving
        ``2 * len(bounds)`` processes. Unlike the unchunked scripts this one
        does not compress its output: the chunks must be merged first, which
        :meth:`get_orbit_library_chunked` does once the script returns.

        Parameters
        ----------
        bounds : list of tuple
            ``(start, number, tag)`` per chunk, from
            :meth:`orbit_chunk_bounds`.

        Returns
        -------
        string
            name of the script written.

        """
        if self.system.is_bar_disk_system():
            orb_prgrm = 'orblib_bar'
        else:
            orb_prgrm = 'orblib_new_mirror'
        # per-process name: concurrently evaluated models sharing this orbit
        # library each write and run their own script, with their own chunk
        # tags in it
        cmd_string = f'cmd_tube_box_orbs_chunked_p{os.getpid()}'
        txt_file = open(cmd_string, 'w')
        txt_file.write('#!/bin/bash\n')
        txt_file.write('# clear flags\n')
        txt_file.write('rm -f datfil/tube_done datfil/box_done '
                       'datfil/tube_box_done\n')
        txt_file.write('# check whether executables exist\n')
        execs = [orb_prgrm]
        if self.LegacyWeightSolver:
            execs += ['triaxmass', 'triaxmass_bar',
                      'triaxmassbin', 'triaxmassbin_bar']
        for f_name in execs:
            txt_file.write(f'test -e {self.legacy_directory}/{f_name} || ' +
                           f'{{ echo "File {self.legacy_directory}/{f_name} ' +
                           'not found." && exit 127; }\n')
        txt_file.write(f'# integrate {len(bounds)} chunk(s) of each orbit '
                       'family, all in parallel\n')
        pids = []
        for fileroot in ('orblib', 'orblibbox'):
            for _, _, tag in bounds:
                stem = f'datfil/{fileroot}{tag}'
                txt_file.write(f'(rm -f {stem}_qgrid.dat.tmp {stem}_qgrid.dat '
                               f'{stem}_losvd_hist.dat\n')
                txt_file.write(f'{self.legacy_directory}/{orb_prgrm} '
                               f'< infil/{fileroot}{tag}.in '
                               f'>> datfil/{fileroot}.log) &\n')
                pid = f'p_{fileroot}{tag}'
                txt_file.write(f'{pid}=$!\n')
                pids.append(pid)
        txt_file.write('# wait for every chunk to finish\n')
        txt_file.write('wait ' + ' '.join(f'${p}' for p in pids) + '\n')
        if self.LegacyWeightSolver:
            # depends on the potential rather than the orbit library, so it is
            # run once here, not once per chunk
            txt_file.write('rm -f datfil/mass_qgrid.dat '
                           'datfil/mass_radmass.dat datfil/mass_aper.dat\n')
            bar = '_bar' if self.system.is_bar_disk_system() else ''
            txt_file.write(f'{self.legacy_directory}/triaxmass{bar} '
                           '< infil/triaxmass.in >> datfil/triaxmass.log\n')
            txt_file.write(f'{self.legacy_directory}/triaxmassbin{bar} '
                           '< infil/triaxmassbin.in '
                           '>> datfil/triaxmassbin.log\n')
        txt_file.close()
        return cmd_string

    def write_executable_for_integrate_orbits_par(self):
        """Write the bash script to calculate orbit libraries
        """
        if self.system.is_bar_disk_system():
            orb_prgrm = 'orblib_bar'
        else:
            orb_prgrm = 'orblib_new_mirror'
        cmd_string = 'cmd_tube_box_orbs'
        txt_file = open(cmd_string, "w")
        txt_file.write('#!/bin/bash\n')
        txt_file.write('# clear flags\n')
        txt_file.write('rm -f datfil/tube_done datfil/box_done '
                       'datfil/tube_box_done\n')
        txt_file.write('# check whether executables exist\n')
        execs = [orb_prgrm]
        if self.LegacyWeightSolver:
            execs += ['triaxmass', 'triaxmass_bar',
                      'triaxmassbin', 'triaxmassbin_bar']
        for f_name in execs:
            txt_file.write(f'test -e {self.legacy_directory}/{f_name} || ' +
                           f'{{ echo "File {self.legacy_directory}/{f_name} ' +
                           'not found." && exit 127; }\n')
        txt_file.write('# now run tube and box orbits in parallel\n')
        txt_file.write('(rm -f datfil/orblib_qgrid.dat.tmp '
                       'datfil/orblib_qgrid.dat '
                       'datfil/orblib_pops.dat datfil/orblib_losvd_hist.dat '
                       'datfil/orblib_pm_hist.dat\n')
        txt_file.write(f'{self.legacy_directory}/{orb_prgrm} < infil/orblib.in '
                        '>> datfil/orblib.log\n')
        if self.LegacyWeightSolver:
            txt_file.write('rm -f datfil/mass_qgrid.dat '
                           'datfil/mass_radmass.dat '
                           'datfil/mass_aper.dat\n')
            bar = '_bar' if self.system.is_bar_disk_system() else ''
            txt_file.write(f'{self.legacy_directory}/triaxmass{bar} '
                '< infil/triaxmass.in >> datfil/triaxmass.log\n')
            txt_file.write(f'{self.legacy_directory}/triaxmassbin{bar} '
                '< infil/triaxmassbin.in >> datfil/triaxmassbin.log\n')
        for f in 'qgrid', 'pm_hist', 'pops', 'losvd_hist':
            f_name = 'datfil/orblib_' + f + '.dat'
            txt_file.write(f'test -e {f_name} '
                           f'&& rm -f {f_name}.bz2 && bzip2 -k {f_name}\n')
            txt_file.write(f'rm -f {f_name}\n')
        txt_file.write('touch datfil/tube_done) &\n')
        txt_file.write('orblib=$!\n')
        txt_file.write('(rm -f datfil/orblibbox_qgrid.dat.tmp '
                       'datfil/orblibbox_qgrid.dat datfil/orblibbox_pops.dat '
                       'datfil/orblibbox_losvd_hist.dat '
                       'datfil/orblibbox_pm_hist.dat\n')
        txt_file.write(f'{self.legacy_directory}/{orb_prgrm} '
                       '< infil/orblibbox.in >> datfil/orblibbox.log\n')
        for f in 'qgrid', 'pm_hist', 'pops', 'losvd_hist':
            f_name = 'datfil/orblibbox_' + f + '.dat'
            txt_file.write(f'test -e {f_name} '
                           f'&& rm -f {f_name}.bz2 && bzip2 -k {f_name}\n')
            txt_file.write(f'rm -f {f_name}\n')
        txt_file.write('touch datfil/box_done) &\n')
        txt_file.write('orblibbox=$!\n')
        txt_file.write('# wait for tube and box orbits to finish\n')
        txt_file.write('wait $orblib $orblibbox\n')
        txt_file.close()
        # returns the name of the executables
        return cmd_string

    def write_executable_for_integrate_orbits(self):
        """Write the bash script to calculate orbit libraries
        """
        if self.system.is_bar_disk_system():
            orb_prgrm = 'orblib_bar'
        else:
            orb_prgrm = 'orblib_new_mirror'
        # tubeorbits
        cmdstr_tube = 'cmd_tube_orbs'
        txt_file = open(cmdstr_tube, "w")
        txt_file.write('#!/bin/bash\n')
        txt_file.write('# clear flags\n')
        txt_file.write('rm -f datfil/tube_done datfil/tube_box_done\n')
        txt_file.write('# check whether executables exist\n')
        execs = [orb_prgrm]
        if self.LegacyWeightSolver:
            execs += ['triaxmass', 'triaxmass_bar',
                      'triaxmassbin', 'triaxmassbin_bar']
        for f_name in execs:
            txt_file.write(f'test -e {self.legacy_directory}/{f_name} || ' +
                           f'{{ echo "File {self.legacy_directory}/{f_name} ' +
                           'not found." && exit 127; }\n')
        txt_file.write('# now run tube orbits\n')
        txt_file.write('rm -f datfil/orblib_qgrid.dat.tmp '
                       'datfil/orblib_qgrid.dat '
                       'datfil/orblib_qgrid.dat.bz2 datfil/orblib_pops.dat '
                       'datfil/orblib_pops.dat.bz2 '
                       'datfil/orblib_losvd_hist.dat '
                       'datfil/orblib_losvd_hist.dat.bz2 '
                       'datfil/orblib_pm_hist.dat '
                       'datfil/orblib_pm_hist.dat.bz2\n')
        txt_file.write(f'{self.legacy_directory}/{orb_prgrm} < infil/orblib.in '
                       '>> datfil/orblib.log\n')
        txt_file.write('rm -f datfil/mass_aper.dat\n')
        if self.LegacyWeightSolver:
            txt_file.write('rm -f datfil/mass_qgrid.dat '
                           'datfil/mass_radmass.dat\n')
            bar = '_bar' if self.system.is_bar_disk_system() else ''
            txt_file.write(f'{self.legacy_directory}/triaxmass{bar} '
                '< infil/triaxmass.in >> datfil/triaxmass.log\n')
            txt_file.write(f'{self.legacy_directory}/triaxmassbin{bar} '
                '< infil/triaxmassbin.in >> datfil/triaxmassbin.log\n')
        for f in 'qgrid', 'pm_hist', 'pops', 'losvd_hist':
            f_name = 'datfil/orblib_' + f + '.dat'
            txt_file.write(f'test -e {f_name} '
                           f'&& bzip2 -kc {f_name} > {f_name}.staging.bz2 '
                           f'&& mv {f_name}.staging.bz2 {f_name}.bz2\n')
            txt_file.write(f'rm -f {f_name}\n')
        txt_file.write('touch datfil/tube_done\n')
        txt_file.close()
        # boxorbits
        cmdstr_box = 'cmd_box_orbs'
        txt_file = open(cmdstr_box, "w")
        txt_file.write('#!/bin/bash\n')
        txt_file.write('# clear flag\n')
        txt_file.write('rm -f datfil/box_done\n')
        txt_file.write('# check whether executable exists\n')
        txt_file.write(f'test -e {self.legacy_directory}/{orb_prgrm} || ' +
                       f'{{ echo "File {self.legacy_directory}/{orb_prgrm} ' +
                       'not found." && exit 127; }\n')
        txt_file.write('# now run box orbits\n')
        txt_file.write('rm -f datfil/orblibbox_qgrid.dat.tmp '
                        'datfil/orblibbox_qgrid.dat '
                        'datfil/orblibbox_qgrid.dat.bz2 '
                        'datfil/orblibbox_pops.dat '
                        'datfil/orblibbox_pops.dat.bz2 '
                        'datfil/orblibbox_losvd_hist.dat '
                        'datfil/orblibbox_losvd_hist.dat.bz2 '
                        'datfil/orblibbox_pm_hist.dat '
                        'datfil/orblibbox_pm_hist.dat.bz2\n')
        txt_file.write(f'{self.legacy_directory}/{orb_prgrm} '
                       '< infil/orblibbox.in >> datfil/orblibbox.log\n')
        for f in 'qgrid', 'pm_hist', 'pops', 'losvd_hist':
            f_name = 'datfil/orblibbox_' + f + '.dat'
            txt_file.write(
                f'test -e {f_name} '
                f'&& bzip2 -kc {f_name} > {f_name}.staging.bz2 '
                f'&& mv {f_name}.staging.bz2 {f_name}.bz2\n')
            txt_file.write(f'rm -f {f_name}\n')
        txt_file.write('touch datfil/box_done\n')
        txt_file.close()
        # returns the name of the executables
        return cmdstr_tube, cmdstr_box

    def read_ics(self):
        # ...
        pass

    def _walk_losvd_records(self, buf, n_pairs):
        """Locate every losvd histogram payload in a raw file buffer.

        Sequential pass over the Fortran record markers of a
        ``{fileroot}_losvd_hist.dat`` file, whose layout after the header
        record is, for each orbit and each aperture::

            record: [int32 ivmin][int32 ivmax]
            record: [float64 x nv]     ONLY present if ivmin <= ivmax

        Parameters
        ----------
        buf : 1d numpy array of uint8
            the whole file, as read by ``np.fromfile``
        n_pairs : int
            number of (orbit, aperture) pairs to read, i.e. norb*n_apertures

        Returns
        -------
        tuple of three 1d int64 arrays, each of length ``n_pairs``:
            start : int32-index of the first float64 of the payload
            ivmin : the record's ivmin
            nv : number of float64 values; ``0`` marks an empty pair
                 (ivmin > ivmax), which has no payload record at all
        """
        n32 = buf.nbytes // 4
        # a memoryview of python ints is ~2x faster to index than a numpy
        # array, and this loop runs once per (orbit, aperture) pair
        mv = memoryview(buf.data)[:n32 * 4].cast('i')
        p = (4 + mv[0] + 4) // 4  # skip the header record
        start = [0] * n_pairs
        ivmin = [0] * n_pairs
        nv = [0] * n_pairs
        for k in range(n_pairs):
            i_min = mv[p + 1]
            i_max = mv[p + 2]
            p += 4
            if i_min <= i_max:
                nbytes = mv[p]
                start[k] = p + 1
                ivmin[k] = i_min
                nv[k] = nbytes >> 3
                p += 2 + (nbytes >> 2)
        return (np.array(start, dtype=np.int64),
                np.array(ivmin, dtype=np.int64),
                np.array(nv, dtype=np.int64))

    @staticmethod
    def _value_batches(nv, max_values):
        """Split pairs into batches holding at most ``max_values`` values each.

        Parameters
        ----------
        nv : 1d numpy array of int
            number of velocity values per (orbit, aperture) pair
        max_values : int
            values per batch. A pair whose own length exceeds this still forms
            a batch by itself, so the split is always progressing.

        Returns
        -------
        1d numpy array of int
            batch boundaries, starting at 0 and ending at ``len(nv)``

        """
        if nv.size == 0:
            return np.array([0])
        cumulative = np.cumsum(nv)
        n_batches = int(cumulative[-1] // max_values) + 1
        # the first pair whose running total passes each multiple of the budget
        cuts = np.searchsorted(
            cumulative, np.arange(1, n_batches) * max_values, side='left') + 1
        return np.unique(np.concatenate(([0], cuts, [nv.size])))

    @staticmethod
    def _gather_float64(buf, byte_offsets):
        """Read float64 values at arbitrary 4-byte-aligned byte offsets.

        Fortran brackets each record with 4-byte length markers, so float64
        payloads are alternately 8-byte aligned and 4-byte misaligned and
        cannot all be read through a single ``view(np.float64)``. Splitting on
        the offset's alignment lets each half be read as an aligned view.

        Parameters
        ----------
        buf : 1d numpy array of uint8
            the file contents
        byte_offsets : 1d numpy array of int
            byte offset of each value to read, each a multiple of 4

        Returns
        -------
        1d numpy array of float64
            the values, in the order their offsets were given

        """
        out = np.empty(byte_offsets.size, dtype=np.float64)
        for phase in (0, 4):
            mask = (byte_offsets & 7) == phase
            if not mask.any():
                continue
            tail = buf[phase:]
            f64 = tail[:tail.size // 8 * 8].view(np.float64)
            out[mask] = f64[(byte_offsets[mask] - phase) >> 3]
        return out

    def _read_losvd_hist_vectorised(self, fname, norb, kin_idx_per_ap,
                                    idx_ap_reset, hist_bins, velhist0,
                                    chunk=2000000):
        """Fill ``velhist0`` from a losvd histogram file, without a read loop.

        Equivalent to the per-record ``FortranFile`` loop in
        :meth:`read_orbit_base`, but reads the file in one go and scatters the
        values with numpy. Results are bit-identical; measured ~24x faster on
        an omega Cen library (12000 orbits x 1012 apertures), where 77% of the
        pairs are empty and the loop spends most of its time reading sentinels.

        Only valid for non-legacy files whose kinematics all have 1d (losvd)
        histograms; :meth:`read_orbit_base` checks this before calling.

        Parameters
        ----------
        fname : string
            the decompressed ``{fileroot}_losvd_hist_{ml}.dat`` file
        norb : int
            number of orbits in the library
        kin_idx_per_ap : 1d numpy array of int
            ``kin_idx_per_ap[i]`` = index of the kinematic set of aperture i
        idx_ap_reset : 1d numpy array of int
            offset subtracted from a global aperture index to get the index
            within its own kinematic set
        hist_bins : list of int
            number of velocity bins, per kinematic set
        velhist0 : list of 3d numpy arrays
            modified in place; element k has shape (norb, nv, n_apertures)
        chunk : int, optional
            number of velocity values processed per batch. The scatter builds
            index arrays as long as the expanded values, so this is what bounds
            peak memory; roughly 56 bytes are needed per value. The default
            keeps the temporaries near 110 MB regardless of library size.
        """
        # memory-map rather than read: the file is only touched once, and this
        # keeps a large library off the heap (and is marginally faster)
        buf = np.memmap(fname, dtype=np.uint8, mode='r')
        n_ap_total = len(kin_idx_per_ap)
        n_pairs = norb * n_ap_total
        start, ivmin, nv = self._walk_losvd_records(buf, n_pairs)
        kin_idx_per_ap = np.asarray(kin_idx_per_ap)
        idx_ap_reset = np.asarray(idx_ap_reset)
        nv0 = np.array([(b - 1) // 2 for b in hist_bins])
        # Split on the number of velocity values, not on the number of pairs:
        # the scatter below builds several arrays as long as the expanded
        # values, so batching by pairs would let peak memory grow with the
        # library's occupancy instead of staying bounded.
        edges = self._value_batches(nv, chunk)
        for lo, hi in zip(edges[:-1], edges[1:]):
            k = np.arange(lo, hi)[nv[lo:hi] > 0]
            if k.size == 0:
                continue
            nvs = nv[k]
            orb = k // n_ap_total
            ap = k - orb * n_ap_total
            kin = kin_idx_per_ap[ap]
            # expand each pair into its nv individual values
            tot = int(nvs.sum())
            within = np.arange(tot) - np.repeat(np.cumsum(nvs) - nvs, nvs)
            vals = self._gather_float64(
                buf, np.repeat(start[k] * 4, nvs) + 8 * within)
            kin_e = np.repeat(kin, nvs)
            # each kinematic set has its own velhist0 array
            for kin_id in np.unique(kin_e):
                m = kin_e == kin_id
                velhist0[kin_id][np.repeat(orb, nvs)[m],
                                 (np.repeat(ivmin[k] + nv0[kin], nvs)
                                  + within)[m],
                                 np.repeat(ap - idx_ap_reset[kin], nvs)[m]] \
                    = vals[m]

    def _read_individual_orbit(self, fort_file, quad_light_grid_sizes):
        """Read individual orbit parameters from file

        Parameters
        ----------
        fort_file : scipy.io.FortranFile object
            The file object to read from, typically pointing to
            {fileroot}.dat (legacy) or {fileroot}_qgrid.dat (new).
        quad_light_grid_sizes : numpy array of shape (4,)
            the orbit library's (l, phi, theta, r) grid sizes

        Returns
        -------
        tuple of length 3:
            orbtypes_dith : 1d numpy array
                the orbit types for each dither in the orbit library
            density_3D_orb : 3d numpy array
                the 3D orbit density
            quad_light : 4d numpy array
                0th to 2nd moments and orbit classification
        """
        _, _, _, _, _ = fort_file.read_ints(np.int32)
        orbtypes_dith = fort_file.read_ints(np.int32)
        quad_light = fort_file.read_reals(float)
        quad_light = np.reshape(quad_light, quad_light_grid_sizes[::-1])
        # quad_light stores orbit features in 3D (r,th,phi) bins.
        # Quad_light[ir,it,ip,XXX] stores
        # - the zeroth moment i.e. density for XXX=0,
        # - the first moments x,y,z,vx,vy,vz for XXX=slice(1,7)
        # - 2nd moments vx^2,vy^2,vz^2,vx*vy,vy*vz,vz*vx for XXX=slice(7,13)
        # - an averaged orbit classification for XXX=slice(13,None)
        # in the bin indexed by (ir,it,ip).
        # We need to extract 3D density for use in weight solving.
        density_3D_orb = quad_light[:,:,:,0]
        return orbtypes_dith, density_3D_orb, quad_light

    def read_orbit_base(self,
                        fileroot,
                        return_intrinsic_moments=False,
                        pops=False):
        """
        Read orbit library from file datfil/{fileroot}.dat.bz2'

        Depending on the DYNAMITE version, the orbit library will reside in
        either datfil/{fileroot}.dat.bz2 (legacy behavior) or in separate
        files datfil/{fileroot}_qgrid.dat.bz2 and one or both of
        datfil/{fileroot}_losvd_hist.dat.bz2 and
        datfil/{fileroot}_pm_hist.dat.bz2 (new behavior).
        If both "legacy" and "new" files exist, default to the new behavior.
        With 'new behavior', populations data (projected masses) may exist in
        datfil/{fileroot}_pops.dat.bz2 and can be read by setting pops=True.

        Parameters
        ----------
        fileroot : string
            This will probably be either 'orblib' or 'orblibbox'.
        return_intrinsic_moments: boolean, optional
            Whether to return_intrinsic_moments of the orblib.
            The default is False.
        pops: boolean, optional
            False: only read the orbit library and velocity
            histograms, as needed by the weight solvers.
            True: only read the population data's orbit densities.
            The default is False.

        Returns
        -------
        If return_intrinsic_moments is False:
        pops==False will return a tuple of type (list, array) where the
        orbit library velocities are stored in the list of Histogram objects,
        and the 3D density of the orbits are stored in the array object.
        pops==True will return a tuple of type (list, None) where the
        populations' projected masses are in the list of Histogram objects.

        return_intrinsic_moments=True will return a tuple (array, list)
        where the array stores the intrinsic moments of the orblib and the
        list contains the bin edges of the 3D grid.

        """
        stars = self.stars
        norb = self.settings['nE'] * self.settings['nI2'] * self.settings['nI3']
        ml = self.parset['ml']
        cur_dir = os.getcwd()
        os.chdir(self.mod_dir)
        check = os.path.isfile(f'datfil/{fileroot}_qgrid.dat.bz2')
        os.chdir(cur_dir)
        legacy_file = False if check else True
        if pops and legacy_file:
            err_msg = f'Pops data not available in legacy mode: {self.mod_dir}.'
            self.logger.error(err_msg)
            raise ValueError(err_msg)
        if pops and return_intrinsic_moments:
            err_msg = 'Pops=True not compatible with return_intrinsic_moments' \
                      f'=True, will set pops to False: {self.mod_dir}.'
            self.logger.warning(err_msg)
            pops = False

        if not pops:  # need orbit properties in 'non-populations' mode only
            os.chdir(self.mod_dir)
            if legacy_file:
                orblib_file = f'datfil/{fileroot}.dat.bz2'
                tmpfname = f'datfil/{fileroot}_{ml}.dat'
            else:
                orblib_file = f'datfil/{fileroot}_qgrid.dat.bz2'
                tmpfname = f'datfil/{fileroot}_qgrid_{ml}.dat'
            # unzip orblib to a temporary file with ml value attached
            # ml value is needed to prevent different processes clashing
            subprocess.run(f'bunzip2 -c {orblib_file} > {tmpfname}', shell=True)
            # read the fortran file
            orblib_in = FortranFile(tmpfname, 'r')
            # read size of orbit library
            # from integrator_setup_write, lines 506 - 5129:
            norb_read, _, _, _, ndith = orblib_in.read_ints(np.int32)
            if norb_read != norb:
                error_msg = f'Number of orbits in {self.mod_dir}{orblib_file}' \
                            f' is {norb_read}, but expected {norb}.'
                self.logger.error(error_msg)
                orblib_in.close()
                os.remove(tmpfname)
                os.chdir(cur_dir)
                raise ValueError(error_msg)
            # from qgrid_setup_write, lines 2339-1350:
            quad_light_grid_sizes = orblib_in.read_ints(np.int32)
            size_ql, size_qph, size_qth, size_qlr = quad_light_grid_sizes
            quad_lr = orblib_in.read_reals(float)
            quad_lth = orblib_in.read_reals(float)
            quad_lph = orblib_in.read_reals(float)
            if return_intrinsic_moments:
                intrinsic_moms = np.zeros((norb,size_qlr,size_qth,size_qph,16))
                intrinsic_grid = [quad_lr, quad_lth, quad_lph]
            orbtypes = np.zeros((norb, ndith**3), dtype=int)
            density_3D = np.zeros((norb, size_qlr, size_qth, size_qph))

            if not legacy_file:  # finish reading orblib_qgrid.dat.bz2
                if return_intrinsic_moments:
                    for j in range(norb):
                        _, _, intrinsic_moms[j] = \
                            self._read_individual_orbit(orblib_in,
                                                        quad_light_grid_sizes)
                else:
                    for j in range(norb):
                        orbtypes[j, :], density_3D[j] , _ = \
                            self._read_individual_orbit(orblib_in,
                                                        quad_light_grid_sizes)
                # done with orblib_qgrid.dat.bz2
                orblib_in.close()
                # remove temporary file
                os.remove(tmpfname)
                os.chdir(cur_dir)
                if return_intrinsic_moments:  # in that case, we are done
                    return intrinsic_moms, intrinsic_grid  ####################
        else:
            density_3D = None

        # next, check if we need to read the losvd_hist and the pm_hist files:
        # note: legacy_file=True never has 2d histograms
        # either pops=False or pops=True and some pops and kins share apertures
        pops_unique = [p for p in stars.population_data if p.kin_aper is None]
        if pops and len(pops_unique) < len(stars.population_data):
            kin_aper_for_pops = \
                [i.kin_aper for i in stars.population_data
                            if i.kin_aper is not None]
        if not pops or len(pops_unique) < len(stars.population_data):
            n_kins = len(stars.kinematic_data)
            hist_widths = [k.hist_width for k in stars.kinematic_data]
            # UNUSED hist_centers = [k.hist_center for k in stars.kinematic_data]
            hist_bins = [k.hist_bins for k in stars.kinematic_data]
            # kinematics k have 1d losvd histograms if type(k.hist_bins)=int and
            # 2d proper motion hists otherwise (k.hist_bins is a 1d array then)
            hist_dim = [1 if type(k.hist_bins) is int else 2
                        for k in stars.kinematic_data]
            if not legacy_file:  # open losvd_hist and pm_hist files if needed
                os.chdir(self.mod_dir)
                if any(i == 1 for i in hist_dim):
                    orblib_file = f'datfil/{fileroot}_losvd_hist.dat.bz2'
                    tmpfname = f'datfil/{fileroot}_losvd_hist_{ml}.dat'
                    subprocess.run(f'bunzip2 -c {orblib_file} > {tmpfname}',
                                shell=True)
                    # read the fortran file
                    orblib_in = FortranFile(tmpfname, 'r')
                if any(i == 2 for i in hist_dim):
                    orblib_file = f'datfil/{fileroot}_pm_hist.dat.bz2'
                    tmpfname_pm = f'datfil/{fileroot}_pm_hist_{ml}.dat'
                    subprocess.run(f'bunzip2 -c {orblib_file} > {tmpfname_pm}',
                                   shell=True)
                    orblib_in_pm = FortranFile(tmpfname_pm, 'r')
            # read the losvd and pm histogram data
            # from histogram_setup_write, lines 1917-1926: (losvd only)
            if legacy_file or any(i == 1 for i in hist_dim):
                _ = orblib_in.read_record(np.int32, np.int32, float)
            # tmp = orblib_in.read_record(np.int32, np.int32, float)
            # nconstr = tmp[0][0] # = total number of apertures for ALL kinematics
            # nvhist = tmp[1][0] # = (nvbins-1)/2 for histo of FIRST kinematic set
            # dvhist = tmp[2][0] # = delta_v in histogram for FIRST kinematic set
            # these nvhist and dvhist are for the first kinematic set only
            # however, orbits are stored N times where N = number of kinematic sets
            # histogram settings for other N-1 sets may be different from the first
            # these aren't stored in orblib.dat so must read from kinematics objects
            self.logger.debug(f'{self.mod_dir}{tmpfname}: '
                              'checking number of velocity bins...')
            hist_bins_flat = [[h] for h in hist_bins if isinstance(h, int)]
            hist_bins_flat += [list(h) for h in hist_bins
                                       if not isinstance(h, int)]
            hist_bins_flat = np.array([i for h in hist_bins_flat for i in h])
            if np.any(hist_bins_flat % 2 == 0):
                error_msg = f'{self.mod_dir}{tmpfname}: all kinematics ' \
                            'need odd number of velocity bins.'
                self.logger.error(error_msg)
                if legacy_file or any(i == 1 for i in hist_dim):
                    orblib_in.close()
                    os.remove(tmpfname)
                if any(i == 2 for i in hist_dim):
                    orblib_in_pm.close()
                    os.remove(tmpfname_pm)
                os.chdir(cur_dir)
                raise ValueError(error_msg)
            self.logger.debug(f'...{self.mod_dir}{tmpfname}: checks ok.')
            n_apertures = [k.n_spatial_bins for k in stars.kinematic_data]
            # get index linking kinematic set to aperture
            # kin_idx_per_ap[i] = N <--> aperture i is from kinematic set N
            kin_idx_per_ap = [np.zeros(n_apertures[i], dtype=int) + i
                              for i in range(n_kins)]
            kin_idx_per_ap = np.concatenate(kin_idx_per_ap)  # is of type int
            # below we loop i_ap from 1-n_total_apertures but will need the
            # index of i_ap for the relevant kinematic set:
            # we use `idx_ap_reset` to do this
            cum_n_apertures = np.cumsum(n_apertures)
            idx_ap_reset = np.concatenate(([0], cum_n_apertures[:-1]))
            # set up a list of arrays to hold the results
            # velhist0 is a list of length n_kins, each element is
            # EITHER a 3D array of shape (norb, nv, na) where nv is the number
            # of 1D velocity bins and na is the number of apertures (=spatial
            # bins) for the kinematic set
            # OR a 4D array of shape (norb, nv[0], nv[1], na) where nv[0] and
            # nv[1] are the number of 2D velocity bins for the proper
            # motion histograms
            velhist0 = []
            for kin_idx, nv in enumerate(hist_bins):
                if hist_dim[kin_idx] == 1:
                    velhist0 += [np.zeros((norb, nv, n_apertures[kin_idx]))]
                elif hist_dim[kin_idx] == 2:
                    velhist0 += [np.zeros((norb,
                                           nv[0],
                                           nv[1],
                                           n_apertures[kin_idx]))]
                else:
                    error_msg = 'Invalid histogram dimension.'
                    self.logger.error(error_msg)    # should never happen
                    raise ValueError(error_msg)
            # Next read the histograms themselves.
            # The loop below does two FortranFile record reads per (orbit,
            # aperture) pair, which for a large library is tens of millions of
            # python-level calls. When the file holds nothing but 1d losvd
            # histograms its layout is regular enough to parse in bulk, so
            # prefer that and keep the loop for everything else.
            vectorised = (not legacy_file
                          and not return_intrinsic_moments
                          and all(i == 1 for i in hist_dim))
            if vectorised:
                self._read_losvd_hist_vectorised(tmpfname,
                                                 norb,
                                                 kin_idx_per_ap,
                                                 idx_ap_reset,
                                                 hist_bins,
                                                 velhist0)
            for j in range(0 if vectorised else norb):
                if legacy_file:  # orbit info is interlaced in the legacy file
                    if return_intrinsic_moments:
                        _, _, intrinsic_moms[j] = \
                            self._read_individual_orbit(orblib_in,
                                                        quad_light_grid_sizes)
                    else:
                        orbtypes[j, :], density_3D[j] , _ = \
                            self._read_individual_orbit(orblib_in,
                                                        quad_light_grid_sizes)
                for i_ap, kin_idx in enumerate(kin_idx_per_ap):
                    i_ap0 = i_ap - idx_ap_reset[kin_idx]
                    if hist_dim[kin_idx] == 1:
                        ivmin, ivmax = orblib_in.read_ints(np.int32)
                        if ivmin <= ivmax:
                            nv0 = (hist_bins[kin_idx] - 1) // 2
                            # ^--- this is an integer since hist_bins is odd
                            tmp = orblib_in.read_reals(float)
                            velhist0[kin_idx][j,
                                              ivmin+nv0:ivmax+nv0+1,
                                              i_ap0] = tmp
                    elif hist_dim[kin_idx] == 2:
                        ivmin0, ivmin1, ivmax0, ivmax1 = \
                            orblib_in_pm.read_ints(np.int32)
                        if ivmin0 <= ivmax0 and ivmin1 <= ivmax1:
                            nv0 = [(hist_bins[kin_idx][0] - 1) // 2,
                                   (hist_bins[kin_idx][1] - 1) // 2]
                            # ^--- these are integers since hist_bins is odd
                            tmp = orblib_in_pm.read_reals(float)
                            tmp = tmp.reshape((ivmax0-ivmin0+1,
                                               ivmax1-ivmin1+1), order='F')
                            velhist0[kin_idx][j,
                                              ivmin0+nv0[0]:ivmax0+nv0[0]+1,
                                              ivmin1+nv0[1]:ivmax1+nv0[1]+1,
                                              i_ap0] = tmp
                    else:
                        error_msg = 'Invalid histogram dimension.'
                        self.logger.error(error_msg)    # should never happen
                        raise ValueError(error_msg)

            if legacy_file or any(i == 1 for i in hist_dim):
                orblib_in.close()
                os.remove(tmpfname)
            if any(i == 2 for i in hist_dim):
                orblib_in_pm.close()
                os.remove(tmpfname_pm)
            os.chdir(cur_dir)
            if return_intrinsic_moments:
                return intrinsic_moms, intrinsic_grid  #######################
            else:
                velhists = []
                for kin_idx, velhist in enumerate(velhist0):
                    if pops and kin_idx not in kin_aper_for_pops:
                        # if pops and this kinematic set is not used by pops
                        continue
                    # if pops is False or this kinematic set is not used by pops
                    if hist_dim[kin_idx] == 1:
                        width = hist_widths[kin_idx]
                        bins = hist_bins[kin_idx]
                        idx_center = (bins-1)/2 # integer since hist_bins is odd
                        idx_center = int(idx_center)
                        dvhist0 = width / bins
                        vedg = np.arange(bins + 1) * dvhist0
                        v = (vedg[1:] + vedg[:-1])/2.
                        v_cent = v[idx_center]
                        vedg -= v_cent
                        velhists += [dyn_kin.Histogram(xedg=vedg, y=velhist)]
                    elif hist_dim[kin_idx] == 2:
                        width = np.array(hist_widths[kin_idx])
                        bins = np.array(hist_bins[kin_idx], dtype=int)
                        # dvhist = width / bins
                        vedg = (np.linspace(start=-width[0] / 2,
                                            stop=width[0] / 2,
                                            num=bins[0] + 1,
                                            endpoint=True),
                                np.linspace(start=-width[1] / 2,
                                            stop=width[1] / 2,
                                            num=bins[1] + 1,
                                            endpoint=True))
                        velhists += [dyn_kin.Histogram2D(xedg=vedg, y=velhist)]
                    else:
                        error_msg = 'Invalid histogram dimension.'
                        self.logger.error(error_msg)    # should never happen
                        raise ValueError(error_msg)
        else:
            velhists = []

        if pops and len(pops_unique) > 0:  # read remaining population data
            n_apertures = [p.n_spatial_bins for p in pops_unique]
            # set up a list of arrays to hold the results
            mass0 = [np.zeros((norb, 1, na)) for na in n_apertures]
            # Next read the projected masses (0d histograms) themselves.
            pops_file = f'datfil/{fileroot}_pops.dat.bz2'
            if not os.path.isfile(self.mod_dir + pops_file):
                error_msg = f'Pops file {self.mod_dir}{pops_file} missing.'
                self.logger.error(error_msg)
                raise FileNotFoundError(error_msg)
            os.chdir(self.mod_dir)
            tmpfname = f'datfil/{fileroot}_pops_{ml}.dat'
            subprocess.run(f'bunzip2 -c {pops_file} > {tmpfname}', shell=True)
            # read the fortran file
            orblib_in = FortranFile(tmpfname, 'r')
            for j in range(norb):
                for mass in mass0:
                    mass[j, 0, :] = orblib_in.read_reals(float)
            orblib_in.close()
            # remove temporary file
            os.remove(tmpfname)
            os.chdir(cur_dir)
            # append the populations data to the velhists
            for mass in mass0:
                vvv = dyn_kin.Histogram(xedg=np.array([-0.5, 0.5]), y=mass)
                velhists += [vvv]
        return velhists, density_3D  #######################

    def duplicate_flip_and_interlace_orblib(self, orblib):
        """flip the tube orbits

        Take an orbit library, create a duplicate library with the velocity
        signs flipped, then interlace the two i.e. so that resulting library
        alternates between flipped/unflipped. This creates an orbit library
        consistent with the Fortran output, enforcing the ordering created by
        the for loops in lines 157-178 of triaxnnls_CRcut.f90

        Parameters
        ----------
        orblib : ``dyn.kinematics.Histogram`` or ``dyn.kinematics.Histogram2D``

        Returns
        -------
        ``dyn.kinematics.Histogram`` or ``dyn.kinematics.Histogram2D``
            the duplicated, flipped and interlaced orblib

        """
        self.logger.debug('Checking for symmetric velocity array...')
        error_msg = 'velocity array must be symmetric'
        if orblib.y.ndim == 3:  # 1D histograms
            assert np.allclose(orblib.xedg, -orblib.xedg[::-1]), error_msg
        else:  # 2D histograms
            assert np.allclose(orblib.xedg[0], -orblib.xedg[0][::-1]), error_msg
            assert np.allclose(orblib.xedg[1], -orblib.xedg[1][::-1]), error_msg
        self.logger.debug('...check ok.')
        vel_d = orblib.y
        if vel_d.ndim == 3:  # 1D histograms (losvd)
            n_orbs, n_vel_bins, n_spatial_bins = vel_d.shape
            reversed_vel_d = vel_d[:, ::-1, :]
            new_vel_d = np.zeros((2*n_orbs, n_vel_bins, n_spatial_bins))
        else:  # 2D histograms (proper motion)
            n_orbs, n_vel_bins1, n_vel_bins2, n_spatial_bins = vel_d.shape
            reversed_vel_d = vel_d[:, ::-1, ::-1, :]
            new_vel_d = np.zeros((2*n_orbs,
                                  n_vel_bins1,
                                  n_vel_bins2,
                                  n_spatial_bins))
        new_vel_d[0::2] = vel_d
        new_vel_d[1::2] = reversed_vel_d
        if vel_d.ndim == 3:  # 1D histograms
            new_orblib = dyn_kin.Histogram(xedg=orblib.xedg, y=new_vel_d)
        else:  # 2D histograms
            new_orblib = dyn_kin.Histogram2D(xedg=orblib.xedg, y=new_vel_d)
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
        """Combine two velocity distribution histograms into one.

        Parameters
        ----------
        orblib1 : ``dyn.kinematics.Histogram`` or ``dyn.kinematics.Histogram2D``
        orblib2 : ``dyn.kinematics.Histogram`` or ``dyn.kinematics.Histogram2D``

        Returns
        -------
        ``dyn.kinematics.Histogram`` or ``dyn.kinematics.Histogram2D``
            the combined orbit libraries

        """
        # check orblibs are compatible
        tmp1 = orblib1.y.shape
        tmp2 = orblib2.y.shape
        n_orbs1 = tmp1[0]
        n_orbs2 = tmp2[0]
        n_spatial_bins1 = tmp1[-1]
        n_spatial_bins2 = tmp2[-1]
        if orblib1.y.ndim == 3:  # 1D histograms
            n_vel_bins1 = tmp1[1]
            n_vel_bins2 = tmp2[1]
        else:  # 2D histograms
            n_vel_bins1 = (tmp1[1:3])
            n_vel_bins2 = (tmp2[1:3])
        self.logger.debug('Checking number of velocity bins...')
        error_msg = 'orblibs have different number of velocity bins'
        assert n_vel_bins1==n_vel_bins2, error_msg
        self.logger.debug('Checking velocity arrays...')
        error_msg = 'orblibs have different velocity arrays'
        if orblib1.y.ndim == 3:  # 1D histograms
            assert np.array_equal(orblib1.x, orblib2.x), error_msg
        else:  # 2D histograms
            assert np.array_equal(orblib1.x[0], orblib2.x[0]), error_msg
            assert np.array_equal(orblib1.x[1], orblib2.x[1]), error_msg
        self.logger.debug('Checking number of spatial bins...')
        error_msg = 'orblibs have different number of spatial bins'
        assert n_spatial_bins1==n_spatial_bins2, error_msg
        self.logger.debug('...checks ok.')
        if orblib1.y.ndim == 3:
            new_vel_d = np.zeros((n_orbs1 + n_orbs2,
                                  n_vel_bins1,
                                  n_spatial_bins1))
        else:
            new_vel_d = np.zeros((n_orbs1 + n_orbs2,
                                  n_vel_bins1[0],
                                  n_vel_bins1[1],
                                  n_spatial_bins1))
        new_vel_d[:n_orbs1] = orblib1.y
        new_vel_d[n_orbs1:] = orblib2.y
        if orblib1.y.ndim == 3:
            new_orblib = dyn_kin.Histogram(xedg=orblib1.xedg, y=new_vel_d)
        else:
            new_orblib = dyn_kin.Histogram2D(xedg=orblib1.xedg, y=new_vel_d)
        return new_orblib

    def read_vel_histograms(self, pops=False):
        """Read the orbit library

        Read box orbits and tube orbits, mirrors the latter, and combines.
        Rescales the velocity axis according to the ``ml`` value. Sets LOSVDs,
        proper motion velocity distributions, and 3D grid/aperture masses of
        the combined orbit library.
        If pops=True, only calculates the populations' projected masses.

        Returns
        -------
        If pops is False, sets the attributes:
            -   ``self.vel_histograms``: a list, whose i'th entry is EITHER a
                ``dyn.kinematics.Histogram`` object holding the orbit library
                LOSVDs OR a ``dyn.kinematics.Histogram2D`` object holding the
                orbit library proper motion velocity distributions, binned for
                the i'th kinematic set
            -   ``self.intrinsic_masses``: 3D grid/intrinsic masses of orbit lib
            -   ``self.projected_masses``: aperture/proj. masses of orbit lib
            -   ``self.n_orbs``: number of orbits in the orbit library
        If pops is True, sets the attribute:
            -   ``self.pops_projected_masses``: aperture/proj. masses of
                populations data

        """
        n_pops = len(self.stars.population_data) if pops else 0

        # TODO: check if this ordering is compatible with weights read in by
        # LegacyWeightSolver.read_weights
        try:
            tube_orblib, tube_density_3D = self.read_orbit_base('orblib',
                                                                pops=pops)
        except:
            self.logger.error('Something went seriously wrong when reading '
                              'the tube orbit library. Check disk quota, file '
                              'integrity, and consistent config files. '
                              f'Model: {self.mod_dir}.')
            raise

        if not self.system.is_bar_disk_system():
            # tube orbits are mirrored/flipped and used twice
            tmp = []
            for t in tube_orblib:
                tmp.append(self.duplicate_flip_and_interlace_orblib(t))
            tube_orblib = tmp
            if n_pops == 0:
                tube_density_3D = np.repeat(tube_density_3D, 2, axis=0)

        # read box orbits
        try:
            box_orblib, box_density_3D = self.read_orbit_base('orblibbox',
                                                              pops=pops)
        except:
            self.logger.error('Something went seriously wrong when reading '
                              'the box orbit library. Check disk quota, file '
                              'integrity, and consistent config files. '
                              f'Model: {self.mod_dir}.')
            raise

        # combine orblibs
        self.vel_histograms = [self.combine_orblibs(t, b)
                               for t, b in zip(tube_orblib, box_orblib)]
        # combine density_3D arrays
        if n_pops == 0:
            density_3D = np.vstack((tube_density_3D, box_density_3D))
            for hist in self.vel_histograms:
                hist.scale_x_values(self.velocity_scaling_factor)
            self.intrinsic_masses = density_3D
            self.n_orbs = self.vel_histograms[0].y.shape[0]
            proj_mass = []
            for hist in self.vel_histograms:
                axes = 1 if hist.y.ndim == 3 else (1, 2)  # 1D or 2D histograms
                proj_mass.append(np.sum(hist.y, axis=axes))
            self.projected_masses = proj_mass
        else:  # note: pops cannot share apertures with 2d histogram kins
            proj_mass = [np.sum(hist.y, axis=1) for hist in self.vel_histograms]
            self.pops_projected_masses = proj_mass

    def read_orbit_intrinsic_moments(self, cache=True):
        """Read the intrinsic moments of the orbit library.

        Moments stored in 3D grid over spherical co-ords (r,theta,phi). This
        function reads the data from files, formats them correctly, and
        converts to physical units.

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
                return_intrinsic_moments=True)
            intmom_tubes = \
                self.duplicate_flip_and_interlace_intmoms(intmom_tubes)
            intmom_boxes, _ = self.read_orbit_base(
                'orblibbox',
                return_intrinsic_moments=True)
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
        ``Vrms^2 = average(vx^2 + vy^2 + vz^2 + 2vx*vy + 2vxvz + 2vyvz)``.
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
        for i in range(norb):
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
        - not a box & ``|Ly|>|Lx|`` & ``|Ly|>|Lz|`` --> y-axis tube "-ish"
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
        # UNUSED idx_box = np.where(bool_box)
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
        # UNUSED n_sum = n_box + n_xtish + n_ytish + n_ztish + n_other
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

    def get_projection_tensor(self,
                              minr=None,
                              maxr=None,
                              r_scale='log',
                              nr=50,
                              nl=61,
                              force_lambda_z=False,
                              dL=1e17):
        projection_tensor_pars = {'minr':minr,
                                  'maxr':maxr,
                                  'r_scale':r_scale,
                                  'nr':nr,
                                  'nl':nl,
                                  'dL':dL,
                                  'force_lambda_z':force_lambda_z}
        self.projection_tensor_pars = projection_tensor_pars
        if r_scale not in ['log', 'linear']:
            txt = f"r_scale must be 'log' or 'linear', not {r_scale}."
            self.logger.error(txt)
            raise ValueError(txt)
        # otherwise, continue...
        if not hasattr(self, 'orb_properties'):
            self.read_orbit_property_file()
        orb_properties = self.orb_properties
        self.classify_orbits(dL=dL, force_lambda_z=force_lambda_z)
        if minr is None:
            minr = np.min(orb_properties['r']).value
        if maxr is None:
            maxr = np.max(orb_properties['r']).value
        log10_r_rng = (np.log10(minr), np.log10(maxr)) if r_scale == 'log' \
                      else (np.nan, np.nan)  # avoid error if minr is zero
        lin_r_rng = (minr, maxr)
        lmd_rng = (-1, 1)
        tot_lmd_rng = (0, 1)
        # store ranges for use in plots
        self.projection_tensor_rng = {
            'log10_r_rng':log10_r_rng,
            'lin_r_rng':lin_r_rng,
            'lmd_rng':lmd_rng,
            'tot_lmd_rng':tot_lmd_rng,
        }
        # store ranges for use in plots
        # log10_r_edg = np.linspace(*log10_r_rng, nr+1)
        # lin_r_edg = np.linspace(*lin_r_rng, nr+1)
        r_edg = np.linspace(*log10_r_rng, nr+1) if r_scale == 'log' \
                else np.linspace(*lin_r_rng, nr+1)
        lmd_edg = np.linspace(*lmd_rng, nl+1)
        tot_lmd_edg = np.linspace(*tot_lmd_rng, nl+1)
        if r_scale == 'log':
            r_idx = np.digitize(np.log10(orb_properties['r'].value),
                                bins=r_edg)
        else:
            r_idx = np.digitize(orb_properties['r'].value, bins=r_edg)
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
