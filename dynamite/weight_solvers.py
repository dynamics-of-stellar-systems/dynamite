import os
import copy
import numpy as np
from astropy import table
from astropy.io import ascii
import subprocess
import logging
from typing import NamedTuple
from scipy import optimize

try:
    import cvxopt
except ModuleNotFoundError:
    pass

try:
    # adelie's rayon thread pool mis-detects core count on high-core-count
    # machines and narrows this process's own CPU affinity to a single core
    # as a side effect of import; every process forked afterwards (pool
    # workers, orblib chunks) would otherwise inherit that one-core mask.
    # Restore what we had rather than opening up to os.cpu_count(): that is
    # the whole machine, so it would also discard a restriction the user
    # meant to impose (taskset is not cgroup-enforced, so nothing else
    # would put it back).
    _saved_affinity = os.sched_getaffinity(0) if hasattr(os, "sched_getaffinity") else None
    import adelie.solver as _adelie_solver

    _ADELIE_AVAILABLE = True
    if _saved_affinity is not None:
        os.sched_setaffinity(0, _saved_affinity)
except ImportError:
    _ADELIE_AVAILABLE = False

from dynamite import constants
from dynamite import analysis
from dynamite import kinematics as dyn_kin


class AdelieProblem(NamedTuple):
    """Everything the adelie path needs, built without materializing A.

    ``X`` has exactly n_rows rows: the sqrt(mu) penalty row REPLACES A's
    total-mass row, which survives only as the small vector ``row0_vec``
    (bitwise equal to A[0], i.e. ones/econ[0])."""

    X: np.ndarray  # (n_rows, n_orbs) F-order, column-scaled
    col_norm: np.ndarray  # (n_orbs,) unit-L2 norms incl. penalty row
    y: np.ndarray  # (n_rows,) ALM target; slot 0 updated per iter
    row0_vec: np.ndarray  # (n_orbs,) == A[0] bitwise
    b0: float  # rhs[0] = total_mass/total_mass_error


def _downcast_orblib(orblib, dtype):
    """Downcast the retained orbit-library data to ``dtype``, in place.

    This is most of the resident memory of a solve, not the matrix itself.
    No-op unless ``dtype`` is float32. ``None`` entries are tolerated so the
    streamed path can call it per set as each one arrives; ``copy=False``
    makes re-converting an already-downcast set free.
    """
    if dtype != np.float32:
        return
    for hist in orblib.vel_histograms:
        if hist is not None:
            hist.y = hist.y.astype(np.float32, copy=False)
    intrinsic = getattr(orblib, "intrinsic_masses", None)
    if intrinsic is not None:
        orblib.intrinsic_masses = intrinsic.astype(np.float32, copy=False)
    projected = getattr(orblib, "projected_masses", None)
    if projected is not None:
        orblib.projected_masses = [
            p.astype(np.float32, copy=False) if p is not None else None for p in projected
        ]


def _scale_columns(X, b_rest, dtype):
    """Unit-L2 scale ``X``'s columns in place; returns ``(col_norm, y)``.

    Shared by both augmented-matrix builders (`_build_augmented_X` from an
    existing A, `construct_adelie_matrix_and_rhs` from the orbit library) so
    they cannot drift. Zero columns are left unscaled.
    """
    col_norm = np.linalg.norm(X, axis=0)
    col_norm[col_norm == 0] = 1.0
    X /= col_norm
    y = np.concatenate([[0.0], b_rest]).astype(dtype)
    return col_norm, y


def chi2_vector_from_residuals(resid_full, row0_sq):
    """Per-row squared residuals of ``A @ w - b`` without materializing A.

    ``resid_full`` holds rows 1..n of that residual (from adelie's
    ``state.resid``, which is exactly the plain residual on rows 1.. because
    ``y[1:] == b_rest``); the total-mass row's contribution arrives separately
    as ``row0_sq = (A[0] @ w - b[0])**2`` because X replaces that row with the
    ALM penalty row. Index 0 of the result is ``row0_sq`` itself, keeping the
    indexing of the A-based chi2_vector (chi2_kin slices
    ``[n_mass_constraints:]``).

    Algebraically identical to ``(A @ w - b)**2``; differs only in rounding,
    since the gemv over A is replaced by the solver's accumulated residual.
    """
    resid_full = np.asarray(resid_full, dtype=np.float64).ravel()
    return np.concatenate(([row0_sq], resid_full[1:] ** 2))


class WeightSolver(object):
    """Generic WeightSolver class

    Specific implementations are defined as sub-classes. Each one should
    have a main method `solve`

    Parameters
    ----------
    config : a ``dyn.config_reader.Configuration`` object
    model : a ``dyn.model.Model`` object
    CRcut : Bool, default False
        whether to use the `CRcut` solution for the counter-rotating orbit
        problem. See Zhu et al. 2018 for more. If `CRcut` is given in the
        configuration file's weight solver settings (which is normally the
        case), this parameter is ignored.

    """

    def __init__(self, config, model, CRcut=False):
        self.logger = logging.getLogger(f"{__name__}.{__class__.__name__}")
        self.config = config
        self.system = config.system
        self.settings = config.settings.weight_solver_settings
        self.model = model
        self.direc_with_ml = model.directory
        self.direc_no_ml = model.directory_noml
        if "CRcut" in self.settings.keys():
            CRcut = self.settings["CRcut"]
        self.CRcut = CRcut
        self.weight_file = f"{self.direc_with_ml}{constants.weight_file}"

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
        weights = 0.0
        chi2_tot = 0.0
        chi2_kin = 0.0
        chi2_kinmap = 0.0
        # ...
        return weights, chi2_tot, chi2_kin, chi2_kinmap

    def chi2_kinmap(self, weights, orblib=None):
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
        orblib : ``dyn.orblib.OrbitLibrary``, optional
            An orbit library whose velocity histograms are UNMUTATED (see the
            note in the code below). Pass one to avoid re-reading the library
            from disk, which dominates the cost of this method. Default is
            None, i.e. read a fresh library.

        Returns
        -------
        chi2_kinmap : float
            chi2 directly calculated from the kinematic maps: sum of
            squared residuals of V, sigma, and GH coefficients from h_3 to h_N

        """
        stars = self.system.get_unique_triaxial_visible_component()
        if any(k.type != "GaussHermite" for k in stars.kinematic_data):
            self.logger.info("All kinematics must be 'GaussHermite' for kinmapchi2. Value set to nan.")
            return float("nan")  # #######################################
        number_gh = self.settings["number_GH"]
        chi2_kinmap = 0.0
        # NOTE: the orbit library used here must have UNMUTATED velocity
        # histograms. construct_nnls_matrix_and_rhs zeroes the first and last
        # velocity bin of each 1D histogram in place, to mimic
        # `triaxnnls_CRcut.f90`; reusing a library still in that state changes
        # chi2_kinmap's value (measured: ~7% off, with weights, chi2_tot and
        # chi2_kin all identical, which is exactly the kind of discrepancy
        # nobody would trace back to here).
        # construct_nnls_matrix_and_rhs restores those two edge slices before
        # returning, so a caller that has already built the NNLS matrix can
        # hand us that same library instead of paying for a second full
        # read - the read is the single most expensive step in a model.
        if orblib is None:
            orblib = self.model.get_orblib()
            orblib.read_vel_histograms()
        for kin_set, kin_data in enumerate(stars.kinematic_data):
            n_gh = min(number_gh, kin_data.max_gh_order)
            coefs = ["v", "sigma"] + [f"h{i}" for i in range(3, n_gh + 1)]
            # get the model's projected masses=flux (unused) and kinematic data
            a = analysis.Analysis(config=self.config, model=self.model, kin_set=kin_set)
            model_gh_coef = a.get_gh_model_kinematic_maps(v_sigma_option="fit", weights=weights, orblib=orblib)
            # get the observed projected masses (unused) and kinematic data
            kinematics_data = kin_data.get_data()
            # calculate chi2_kinmap
            for coef in coefs:
                obs_val = np.array(kinematics_data[coef])
                mod_val = np.array(model_gh_coef[coef])
                err_val = np.array(kinematics_data["d" + coef])
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
        self.logger = logging.getLogger(f"{__name__}.{__class__.__name__}")
        self.legacy_directory = self.config.settings.legacy_settings["directory"]
        self.sformat = self.system.parameters[0].sformat  # this is ml's format
        ml_idx = self.direc_with_ml.rindex("/ml")
        ml_str = self.direc_with_ml[ml_idx + 3 :]
        self.ml = float(ml_str[:-1]) if ml_str[-1] == "/" else float(ml_str)
        self.fname_nn_kinem = self.direc_with_ml + "nn_kinem.out"
        self.fname_nn_nnls = self.direc_with_ml + "nn_nnls.out"
        # check the format of the orbit library files
        # check == True means there are 2 orblib_* and 2 orblibbox_* files
        # check == False means there are only two orblib files,
        # orblib.dat.bz2 and orbibbox.dat.bz2 (legacy behavior)
        pth = self.direc_no_ml + "datfil/"
        check = (
            os.path.isfile(f"{pth}orblib_qgrid.dat.bz2")
            and os.path.isfile(f"{pth}orblib_losvd_hist.dat.bz2")
            and os.path.isfile(f"{pth}orblibbox_qgrid.dat.bz2")
            and os.path.isfile(f"{pth}orblibbox_losvd_hist.dat.bz2")
        )
        self.legacy_files = False if check else True
        # prepare fortran input file for nnls
        self.copy_kinematic_data()
        self.create_fortran_input_nnls()
        self.logger.info(
            f"{__class__.__name__} is DEPRECATED and "
            "will be removed in a future version of "
            f"DYNAMITE. Use weight solver type NNLS instead "
            f"of {__class__.__name__} if you can."
        )

    def copy_kinematic_data(self):
        """Copy kin data to infil/ direc"""
        stars = self.system.get_unique_triaxial_visible_component()
        kinematics = stars.kinematic_data
        # convert kinematics to old format to input to fortran
        for i in np.arange(len(kinematics)):
            if len(kinematics) == 1:
                old_filename = self.direc_no_ml + "infil/kin_data.dat"
            else:
                old_filename = self.direc_no_ml + "infil/kin_data_" + str(i) + ".dat"
            kinematics[i].convert_to_old_format(old_filename, self.settings)
        # combine all kinematics into one file
        if len(kinematics) > 1:
            if not all(isinstance(kin, dyn_kin.GaussHermite) for kin in kinematics):
                text = "Multiple kinematics: all must be GaussHermite"
                self.logger.error(text)
                raise ValueError(text)
            # make a dummy 'kins_combined' object ...
            kins_combined = copy.deepcopy(kinematics[0])
            # ...replace data attribute with stacked table of all kinematics
            kins_combined.data = table.vstack([k.get_data() for k in kinematics])
            kins_combined.n_apertures = len(kins_combined.data)
            kins_combined.max_gh_order = self.settings["number_GH"]
            old_filename = self.direc_no_ml + "infil/kin_data_combined.dat"
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
        ml_scaling_factor = self.config.all_models.get_model_velocity_scaling_factor(model=self.model)
        # -------------------
        # write nn.in
        # -------------------
        n_kin = len(self.system.get_unique_triaxial_visible_component().kinematic_data)

        if n_kin == 1:
            kin_data_file = "kin_data.dat"

        else:
            kin_data_file = "kin_data_combined.dat"

        text = (
            "infil/parameters_pot.in"
            + "\n"
            + str(self.settings["regularisation"])
            + "                                  [ regularization strength, 0 = no regularization ]"
            + "\n"
            + f"ml{self.ml:{self.sformat}}/nn\n"
            + "datfil/mass_qgrid.dat"
            + "\n"
            + "datfil/mass_aper.dat"
            + "\n"
            + str(self.settings["number_GH"])
            + "	                           [ # of GH moments to constrain the model]"
            + "\n"
            + "infil/"
            + kin_data_file
            + "\n"
            + str(self.settings["lum_intr_rel_err"])
            + "                               [ relative error for intrinsic luminosity ]"
            + "\n"
            + str(self.settings["sb_proj_rel_err"])
            + "                               [ relative error for projected SB ]"
            + "\n"
            + str(ml_scaling_factor)
            + "                                [ scale factor related to M/L, sqrt( (M/L)_k / (M/L)_ref ) ]"
            + "\n"
        )
        if self.legacy_files:
            text += 2 * f"datfil/orblib_{self.ml}.dat\n" + 2 * f"datfil/orblibbox_{self.ml}.dat\n"  # yes, really...
        else:
            for f in "_qgrid", "_losvd_hist", "box_qgrid", "box_losvd_hist":
                text += f"datfil/orblib{f}_{self.ml}.dat\n"
        text += str(self.settings["nnls_solver"]) + "                                  [ nnls solver ]"

        nn_file = open(self.direc_no_ml + f"ml{self.ml:{self.sformat}}/nn.in", "w")
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
            self.logger.info(f"Reading NNLS solution from existing output {self.weight_file}.")
            results = ascii.read(self.weight_file)
            weights = results["weights"]
            chi2_tot = results.meta["chi2_tot"]
            chi2_kin = results.meta["chi2_kin"]
            chi2_kinmap = results.meta["chi2_kinmap"]
        else:
            fname_nn_orbmat = self.direc_with_ml + "nn_orbmat.out"
            # If legacy result files do not exist, run weight solving.
            check = (
                os.path.isfile(self.fname_nn_kinem)
                and os.path.isfile(self.fname_nn_nnls)
                and os.path.isfile(fname_nn_orbmat)
            )
            if ignore_existing_weights or not check:
                for f in [self.fname_nn_kinem, self.fname_nn_nnls, fname_nn_orbmat]:
                    if os.path.isfile(f):
                        os.remove(f)
                # set the current directory to the directory in which
                # the models are computed
                cur_dir = os.getcwd()
                os.chdir(self.direc_no_ml)
                cmdstr = self.write_executable_for_weight_solver()
                with open(cmdstr) as f:
                    for line in f:
                        i = line.find(">>")
                        if i >= 0:
                            j = line.find(".log")
                            logfile = line[i + 3 : j + 4]
                            break
                self.logger.info("Fitting orbit library to the kinematic " + f"data: {logfile[: logfile.rindex('/')]}")
                p = subprocess.run(
                    "bash " + cmdstr,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    shell=True,
                )
                # clean up decompressed files
                for f_name in [
                    f"datfil/orblib_{self.ml}.dat",
                    f"datfil/orblibbox_{self.ml}.dat",
                ]:
                    if os.path.isfile(f_name):
                        os.remove(f_name)
                log_file = f"Logfile: {self.direc_no_ml + logfile}."
                if not p.stdout.decode("UTF-8"):
                    self.logger.info(f"...done, NNLS problem solved - {cmdstr} exit code {p.returncode}. {log_file}")
                else:
                    text = f"...failed! {cmdstr} exit code {p.returncode}. Message: {p.stdout.decode('UTF-8')}"
                    if p.returncode == 127:  # command not found
                        text += "Check DYNAMITE legacy_fortran executables."
                        self.logger.error(text)
                        os.chdir(cur_dir)
                        raise FileNotFoundError(text)
                    text += f"{log_file} Be wary: DYNAMITE may crash..."
                    self.logger.warning(text)
                    os.chdir(cur_dir)
                    raise RuntimeError(text)
                # set the current directory to the dynamite directory
                os.chdir(cur_dir)
                # delete existing .yaml files and copy current config file
                # into model directory
                self.config.copy_config_file(self.direc_with_ml)
            else:  # If legacy output files exist, just create the weight file
                self.logger.info("Reading NNLS solution from existing legacy output and converting to weights file.")
            # Now the legacy result files exist -> read, calculate
            # kinmapchi2, and save to the weight file.
            weights, chi2_tot, chi2_kin = self.get_weights_and_chi2_from_orbmat_file()
            chi2_kinmap = self.chi2_kinmap(weights)
            # save the output
            results = table.Table()
            results["weights"] = weights
            results.meta = {
                "chi2_tot": chi2_tot,
                "chi2_kin": chi2_kin,
                "chi2_kinmap": chi2_kinmap,
            }
            results.write(self.weight_file, format="ascii.ecsv", overwrite=True)
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
        nn = f"ml{self.ml:{self.sformat}}/nn"
        cmdstr = f"cmd_nnls_{self.ml}"
        txt_file = open(cmdstr, "w")
        txt_file.write("#!/bin/bash" + "\n")
        txt_file.write("# if the gzipped orbit library exist unzip it" + "\n")
        txt_file.write(
            f"test -e datfil/orblib_{self.ml}.dat || bunzip2 -c  datfil/orblib.dat.bz2 > datfil/orblib_{self.ml}.dat"
            + "\n"
        )
        txt_file.write(
            f"test -e datfil/orblibbox_{self.ml}.dat || bunzip2 -c  datfil/orblibbox.dat.bz2 > datfil/orblibbox_{self.ml}.dat"
            + "\n"
        )
        if self.system.is_bar_disk_system():
            txt_file.write(
                f"test -e {self.legacy_directory}/triaxnnls_bar"
                + f' || {{ echo "File {self.legacy_directory}/triaxnnls_bar not found." && exit 127; }}\n'
            )
            txt_file.write(
                "test -e "
                + str(nn)
                + "_kinem.out || "
                + self.legacy_directory
                + f"/triaxnnls_bar < {nn}.in >> {nn}ls.log "
                "|| exit 1\n"
            )
        elif self.CRcut is True:
            txt_file.write(
                f"test -e {self.legacy_directory}/triaxnnls_CRcut"
                + f' || {{ echo "File {self.legacy_directory}/triaxnnls_CRcut not found." && exit 127; }}\n'
            )
            txt_file.write(
                "test -e "
                + str(nn)
                + "_kinem.out || "
                + self.legacy_directory
                + f"/triaxnnls_CRcut < {nn}.in >> {nn}ls.log "
                "|| exit 1\n"
            )
        else:
            txt_file.write(
                f"test -e {self.legacy_directory}/triaxnnls_noCRcut"
                + f' || {{ echo "File {self.legacy_directory}/triaxnnls_noCRcut not found." && exit 127; }}\n'
            )
            txt_file.write(
                "test -e "
                + str(nn)
                + "_kinem.out || "
                + self.legacy_directory
                + f"/triaxnnls_noCRcut < {nn}.in >> {nn}ls.log "
                "|| exit 1\n"
            )
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
        fname = self.direc_with_ml + "nn_orb.out"
        col_names = [
            "orb_idx",
            "E_idx",
            "I2_idx",
            "I3_idx",
            "totalnotregularizable",  # see line 535 of orblib_f.f90
            "orb_type",
            "weight",
            "lcut",
        ]  # lines 1321-1322 of triaxnnls_CRcut.f90
        # NOTE: column 'lcut' is not present if different "triaxnnls" file used
        dtype = [int, int, int, int, int, int, float, int]
        weights = np.genfromtxt(fname, skip_header=1, names=col_names, dtype=dtype)
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
        fname = self.direc_with_ml + "nn_orbmat.out"
        orbmat_shape = np.loadtxt(fname, max_rows=1, dtype=int)
        orbmat_size = np.prod(orbmat_shape)
        tmp = np.loadtxt(fname, skiprows=1)
        orbmat = tmp[0:orbmat_size]
        orbmat = np.reshape(orbmat, orbmat_shape)
        orbmat = orbmat.T
        rhs = tmp[orbmat_size : orbmat_size + orbmat_shape[1]]
        solution = tmp[orbmat_size + orbmat_shape[1] :]
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
        chi2_vector = (np.dot(A, weights) - b) ** 2.0
        chi2_tot = np.sum(chi2_vector)

        stars = self.system.get_unique_triaxial_visible_component()
        if self.system.is_bar_disk_system():
            mge = stars.mge_lum + stars.disk_lum
        else:
            mge = stars.mge_lum

        intrinsic_masses = mge.get_intrinsic_masses_from_file(self.direc_no_ml)
        projected_masses = mge.get_projected_masses_from_file(self.direc_no_ml)
        n_intrinsic = np.prod(intrinsic_masses.shape)
        n_apertures = len(projected_masses)
        chi2_kin = np.sum(chi2_vector[1 + n_intrinsic + n_apertures :])
        return weights, chi2_tot, chi2_kin

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
        # nvel = np.int64(a[0])
        # ncon = np.int64(a[0])
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
        lines = [line.rstrip("\n").split() for line in open(infile)]
        output = []
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
        one of ``scipy``, ``cvxopt`` or ``adelie``. ``adelie`` uses
        coordinate-descent BVLS with an augmented Lagrangian on the total-mass
        constraint; see ``solve_adelie_alm``

    Weight solver settings
    ----------------------
    nnls_dtype : string, optional
        ``'float64'`` (default) or ``'float32'``. ``'float32'`` roughly
        halves the memory of the retained orbit-library data
        (``vel_histograms``/``intrinsic_masses``/``projected_masses``) and
        of the NNLS matrix/solve arrays, at the cost of a small precision
        loss. Validated against float64 on NGC6278 (matrix-only and
        end-to-end): chi2 agrees to <0.001%, KKT violation stays within the
        range of previously accepted float64 solutions. Not yet validated on
        datasets whose constraint-row scaling differs from NGC6278/omega
        Cen - carries the same caveat as ``adelie_mu`` (see
        ``solve_adelie_alm``). Applies to the orbit-library data and the
        constructed matrix regardless of ``nnls_solver``, but has only been
        validated with ``nnls_solver='adelie'``.

    """

    def __init__(self, nnls_solver=None, **kwargs):
        super().__init__(**kwargs)
        self.logger = logging.getLogger(f"{__name__}.{__class__.__name__}")
        if nnls_solver is None:
            nnls_solver = self.settings["nnls_solver"]
        assert nnls_solver in ["scipy", "cvxopt", "adelie"], "Unknown nnls_solver"
        self.nnls_solver = nnls_solver
        # ALM settings for the adelie solver. On NGC6278, where scipy provides
        # a verified optimum, mu between 1e5 and 1e7 reproduces it; 1e7 gave the
        # lowest KKT violation (7e-11) and a monotone gap, whereas 1e6
        # oscillated. See docs/source/adelie_branch_migration.md.
        self.adelie_mu = float(self.settings.get("adelie_mu", 1.0e7))
        self.adelie_alm_iters = int(self.settings.get("adelie_alm_iters", 200))
        self.adelie_tol = float(self.settings.get("adelie_tol", 1.0e-10))
        self.adelie_gap_tol = float(self.settings.get("adelie_gap_tol", 1e-10))
        # Optional float32 mode: roughly halves the memory of the orbit
        # library data retained for the solve (vel_histograms/intrinsic_masses
        # /projected_masses) and of the NNLS matrix/solve arrays. Validated
        # against float64 on NGC6278 (matrix-only and end-to-end): chi2 agrees
        # to <0.001%, KKT violation stays within the range of previously
        # accepted float64 solutions. Not yet validated on datasets whose
        # constraint-row scaling differs from NGC6278/omega Cen - adelie_mu
        # itself carries the same caveat (see solve_adelie_alm docstring).
        nnls_dtype = self.settings.get("nnls_dtype", "float64")
        assert nnls_dtype in ("float32", "float64"), "nnls_dtype must be 'float32' or 'float64'"
        self.nnls_dtype = np.float32 if nnls_dtype == "float32" else np.float64
        # Stream orbit-library histogram reads set-by-set in the fused adelie
        # constructor. Pure memory setting: results are bit-identical either
        # way (validated by dev_tests/_real_fused_check.py).
        self.stream_reads = bool(self.settings.get("stream_orblib_reads", False))
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
                - ``self.total_mass``
                - ``self.total_mass_error``
                -   constraint counts ``self.n_intrinsic``, ``self.n_apertures``
                    and ``self.n_mass_constraints``

        """
        if self.system.is_bar_disk_system():
            mge = self.system.get_unique_bar_component().mge_lum_tot
        else:
            mge = self.system.get_unique_triaxial_visible_component().mge_lum
        # intrinsic mass
        self.intrinsic_masses = mge.get_intrinsic_masses(self.model, nocalc=True)[1]
        self.intrinsic_mass_error = self.settings["lum_intr_rel_err"]
        # projected
        self.projected_masses = mge.get_projected_masses(nocalc=True)
        self.projected_mass_error = self.settings["sb_proj_rel_err"]
        # total mass constraint
        self.total_mass = np.sum(self.intrinsic_masses)
        self.total_mass_error = max(abs(1.0 - self.total_mass), 1e-8)
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
            (orbmat, rhs). ``orbmat`` has shape (n_constraints, n_orbs) and is
            **Fortran-ordered** - see the note below, and keep it that way.

        Notes
        -----
        This assembles the largest array in the run (125 GiB for omega Cen at
        45000 orbits), so it is written to touch that array as few times as
        possible: sizes are computed up front rather than grown by np.vstack,
        the memory order is chosen to match how the blocks arrive and how the
        solvers want to read them, and blocks are written through a reshaped
        view of the destination so no full-size temporary is ever materialised.
        Each of those is commented at its site.

        On the NGC5139 MUSE+HST test library (322403 x 3840, 9.9 GiB) this
        takes 3.9s, against 17.2s before those three changes, and produces a
        bit-identical matrix. dev_tests/_real_orblib_check.py reproduces both
        the timing and the comparison; dev_tests/test_nnls_matrix_assembly.py
        checks the assembly arithmetic without needing a stored orbit library.

        """
        # construct vector of observed constraints (con), errors (econ) and
        # matrix of orbit propertites (orbmat)
        dtype = self.nnls_dtype
        stars = self.system.get_unique_triaxial_visible_component()
        # Size con/econ/orbmat up front. The observed values depend only on the
        # kinematic data, not on the orbit library, so this pre-pass is cheap -
        # and it avoids growing orbmat by np.vstack per kinematic set, which
        # reallocates and copies the whole matrix each time (~125 GiB for
        # omega Cen). The results are kept and reused in the loop below.
        obs_values = [kins.get_observed_values_and_uncertainties(self.settings) for kins in stars.kinematic_data]
        n_rows = self.n_mass_constraints + sum(np.size(v) for v, _ in obs_values)
        con = np.zeros(n_rows, dtype=dtype)
        econ = np.zeros(n_rows, dtype=dtype)
        # F order matters here. Each kinematic block arrives as
        # (n_orbs, n_constraints) and is written in transposed; against an
        # F-ordered destination that is a memcpy, against a C-ordered one it is
        # a strided shuffle of the whole matrix (measured 1.00s -> 0.08s per
        # 2.3 GiB block). It also hands the solvers the layout they want:
        # solve_adelie_alm's X is F-contiguous, so building it from A stops
        # being a full reorder (2.38s -> 0.33s). Only the chi2 matvec A @ w is
        # slightly slower (0.04s -> 0.07s), which is noise beside the rest.
        orbmat = np.zeros((n_rows, orblib.n_orbs), dtype=dtype, order="F")
        # total mass
        con[0] = self.total_mass
        econ[0] = self.total_mass_error
        orbmat[0, :] = 1.0
        # intrinsic mass
        idx = slice(1, 1 + self.n_intrinsic)
        con[idx] = np.ravel(self.intrinsic_masses)
        error = self.intrinsic_masses * self.intrinsic_mass_error
        error = np.abs(np.ravel(error))
        error[np.where(error <= 0.0)] = 1.0e-16
        econ[idx] = np.abs(np.ravel(error))
        orb_int_masses = orblib.intrinsic_masses
        orb_int_masses = np.reshape(orb_int_masses, (orblib.n_orbs, -1))
        orbmat[idx, :] = orb_int_masses.T
        # projected mass
        idx = slice(1 + self.n_intrinsic, 1 + self.n_intrinsic + self.n_apertures)
        con[idx] = self.projected_masses
        econ[idx] = np.abs(self.projected_masses * self.projected_mass_error)
        orbmat[idx, :] = np.hstack(orblib.projected_masses).T
        # add kinematics to con, econ, orbmat
        kins_and_orb_veldist = zip(stars.kinematic_data, orblib.vel_histograms, obs_values)
        idx_ap_start = 0
        idx_row = self.n_mass_constraints
        for kins, orb_veldist, tmp in kins_and_orb_veldist:
            # pick out the projected masses for this kinematic set
            n_ap = kins.n_spatial_bins  # OK for all kinematics
            idx_ap_end = idx_ap_start + n_ap
            prj_mass_i = self.projected_masses[idx_ap_start:idx_ap_end]
            idx_ap_start += n_ap
            obs_kins, obs_kins_err, orb_kins = self._prepare_kinematic_block(kins, orb_veldist, tmp, prj_mass_i)
            # append constraints/errors/orbits to con/econ/orbmat
            n_orb_constraints = orb_kins.size // orblib.n_orbs
            idx_row_end = idx_row + obs_kins.size
            assert n_orb_constraints == obs_kins.size, (
                f"{type(kins).__name__}: orbit library gives "
                f"{n_orb_constraints} constraints per orbit but the "
                f"kinematic data gives {obs_kins.size}"
            )
            # slice assignment casts to dtype in place, no extra full copy
            con[idx_row:idx_row_end] = obs_kins
            econ[idx_row:idx_row_end] = np.ravel(obs_kins_err)
            # orb_kins comes back from transform_orblib_to_observables as a
            # transposed/moveaxis'd view, so np.reshape(orb_kins, (n_orbs, -1))
            # cannot be a view and copies the whole block - 2.1s and a full
            # extra copy of the matrix on NGC5139, so ~100 GiB for omega Cen.
            # Reshaping the *destination* is a view instead, and the write then
            # goes straight from orb_kins with no temporary at all.
            #
            # ndarray.reshape() silently returns a COPY when it cannot return a
            # view. As an assignment target that is a silent no-op: the block
            # would be left as zeros with no error raised anywhere. The assert
            # is what makes that failure loud, so do not drop it. (The .shape=
            # setter raises natively, but is deprecated in numpy 2.5 and we
            # support numpy>=1.26.)
            dest = orbmat[idx_row:idx_row_end, :].T.reshape(orb_kins.shape)
            assert np.shares_memory(dest, orbmat), (
                "orbmat block write got a copy, not a view - block would be silently left at zero"
            )
            dest[...] = orb_kins
            idx_row = idx_row_end
        # divide constraint vector and matrix by errors
        if np.any(con[econ == 0] != 0):
            txt = "Weight solving fail: zero errors for nonzero constraints!"
            self.logger.error(txt)
            raise ValueError(txt)
        # previous statement: rhs = con/econ, np.divide has the "where" clause
        rhs = np.zeros_like(con)
        np.divide(con, econ, out=rhs, where=econ != 0)  # con = econ = 0 is ok
        if np.any(np.ravel(orbmat[econ == 0]) != 0):
            err_loc = np.nonzero(((orbmat != 0).T * (econ == 0)).T)
            txt = (
                f"Weight solving problem in {self.direc_with_ml}: "
                "zero errors for nonzero matrix coefficients at "
                f"[constraint no, orbit no] = {err_loc}! Matrix value(s) "
                f"there ({orbmat[err_loc]}) will be considered zero."
            )
            self.logger.warning(txt)
            orbmat[err_loc] = 0
        # previous statement: orbmat = (orbmat.T/econ).T, np.divide has "where"
        orbmat = orbmat.T
        np.divide(orbmat, econ, out=orbmat, where=econ != 0)
        return orbmat.T, rhs

    def construct_adelie_matrix_and_rhs(self, orblib):
        """Assemble adelie's augmented matrix directly — A never exists.

        Writes the sqrt(mu) penalty row into X row 0, streams every
        constraint block straight into rows 1.., divides by econ, then
        finishes with the SAME col_norm/divide/y steps as
        _build_augmented_X, so X/col_norm/y are bit-identical to the
        two-step construction used by the scipy/cvxopt paths. Saves one full
        matrix of RAM (~125 GiB for omega Cen), which is what lets several
        weight solves share a node.

        Returns an :class:`AdelieProblem`.
        """
        dtype = self.nnls_dtype
        stars = self.system.get_unique_triaxial_visible_component()
        # observed values depend only on the kinematic data; same pre-pass,
        # kept and reused, as in construct_nnls_matrix_and_rhs
        obs_values = [kins.get_observed_values_and_uncertainties(self.settings) for kins in stars.kinematic_data]
        n_rows = self.n_mass_constraints + sum(np.size(v) for v, _ in obs_values)
        sqrt_mu = np.sqrt(self.adelie_mu)
        con = np.zeros(n_rows, dtype=dtype)
        econ = np.zeros(n_rows, dtype=dtype)

        if self.stream_reads:
            # Stream one kinematic set at a time: read -> transform -> write
            # its block into X -> free, so only that set's histograms and the
            # accumulated X are ever co-resident.
            orblib.read_vel_histograms(kin_sets=[0], skip_density=False)
            _downcast_orblib(orblib, self.nnls_dtype)
        n_orbs = orblib.n_orbs
        X = np.empty((n_rows, n_orbs), dtype=dtype, order="F")
        # row 0 is the ALM penalty row; it REPLACES A's total-mass row
        X[0, :] = sqrt_mu
        con[0] = self.total_mass
        econ[0] = self.total_mass_error
        # intrinsic mass -> X rows 1..n_intrinsic
        idx = slice(1, 1 + self.n_intrinsic)
        con[idx] = np.ravel(self.intrinsic_masses)
        error = np.abs(np.ravel(self.intrinsic_masses * self.intrinsic_mass_error))
        error[np.where(error <= 0.0)] = 1.0e-16
        econ[idx] = error
        X[idx, :] = np.reshape(orblib.intrinsic_masses, (n_orbs, -1)).T
        # projected mass: the constraint vector is known up front, the matrix
        # rows are filled per set inside the loop below (each set owns a
        # contiguous aperture range, so this is the same content the old
        # np.hstack(projected_masses).T produced, without the full-size
        # temporary it needed to build it).
        idx_prj = slice(1 + self.n_intrinsic, 1 + self.n_intrinsic + self.n_apertures)
        con[idx_prj] = self.projected_masses
        econ[idx_prj] = np.abs(self.projected_masses * self.projected_mass_error)
        # kinematics: identical block sequence in both read modes. Rows are
        # independent, so streaming produces a bit-identical X.
        idx_ap_start = 0
        idx_row = self.n_mass_constraints
        for si, kins in enumerate(stars.kinematic_data):
            if self.stream_reads and si > 0:
                orblib.read_vel_histograms(kin_sets=[si], skip_density=True)
                _downcast_orblib(orblib, self.nnls_dtype)
            orb_veldist = orblib.vel_histograms[si]
            assert orb_veldist is not None, f"no histogram for kinematic set {si}"
            n_ap = kins.n_spatial_bins
            prj_mass_i = self.projected_masses[idx_ap_start : idx_ap_start + n_ap]
            prj_parts_i = orblib.projected_masses[si]
            assert prj_parts_i is not None, f"no projected masses for kinematic set {si}"
            X[1 + self.n_intrinsic + idx_ap_start : 1 + self.n_intrinsic + idx_ap_start + n_ap, :] = prj_parts_i.T
            obs_kins, obs_kins_err, orb_kins = self._prepare_kinematic_block(
                kins, orb_veldist, obs_values[si], prj_mass_i
            )
            n_orb_constraints = orb_kins.size // n_orbs
            idx_row_end = idx_row + obs_kins.size
            assert n_orb_constraints == obs_kins.size, (
                f"{type(kins).__name__}: orbit library gives "
                f"{n_orb_constraints} constraints per orbit but the "
                f"kinematic data gives {obs_kins.size}"
            )
            con[idx_row:idx_row_end] = obs_kins
            econ[idx_row:idx_row_end] = obs_kins_err
            dest = X[idx_row:idx_row_end, :].T.reshape(orb_kins.shape)
            assert np.shares_memory(dest, X), (
                "X block write got a copy, not a view - block would be silently left at zero"
            )
            dest[...] = orb_kins
            idx_row = idx_row_end
            idx_ap_start += n_ap
            if self.stream_reads:
                # free this set's histograms before reading the next one;
                # glibc returns these mmap-backed pages immediately
                orblib.vel_histograms[si] = None
        # guards: mirror the stock constructor. Row 0 has no econ semantics
        # (it is the penalty row); A's total-mass row becomes row0_vec below.
        if np.any(con[econ == 0] != 0):
            txt = "Weight solving fail: zero errors for nonzero constraints!"
            self.logger.error(txt)
            raise ValueError(txt)
        rhs = np.zeros_like(con)
        np.divide(con, econ, out=rhs, where=econ != 0)  # con = econ = 0 ok
        econ_body = econ[1:]
        # Only zero-error rows can offend, and normally there are none. Restrict
        # to those rows instead of building a (n_rows-1, n_orbs) bool mask: at
        # production scale that mask is ~16 GiB, and the `&` allocates a second
        # one before the first is freed - ~31 GiB transient at the exact moment
        # assembly is already at its RSS peak. The stock constructor above gets
        # this for free by masking rows first; keep the fused path equivalent.
        bad = np.nonzero(econ_body == 0)[0]
        if bad.size:
            blk = X[1:, :][bad]  # (n_bad, n_orbs), n_bad is normally 0
            nz_rows, nz_cols = np.nonzero(blk)
            if nz_rows.size:
                rr = bad[nz_rows]
                txt = (
                    f"Weight solving problem in {self.direc_with_ml}: "
                    "zero errors for nonzero matrix coefficients at "
                    f"[constraint no, orbit no] = {(rr + 1, nz_cols)}! Matrix "
                    f"value(s) there ({blk[nz_rows, nz_cols]}) will be "
                    "considered zero."
                )
                self.logger.warning(txt)
                X[1 + rr, nz_cols] = 0
        # divide rows by their errors: the same elementwise op as the stock
        # transposed-view divide, restricted to rows 1.. (row 0 of A becomes
        # row0_vec below, divided elementwise by econ[0] exactly as stock's
        # broadcast divide did to it).
        body = X[1:, :].T  # view (n_orbs, n_rows-1)
        np.divide(body, econ_body, out=body, where=econ_body != 0)
        col_norm, y = _scale_columns(X, rhs[1:], dtype)
        row0_vec = np.full(n_orbs, 1.0, dtype=dtype) / econ[0]
        return AdelieProblem(X=X, col_norm=col_norm, y=y, row0_vec=row0_vec, b0=float(rhs[0]))

    def _prepare_kinematic_block(self, kins, orb_veldist, tmp, prj_mass_i):
        """Shared by all NNLS constructors.

        Scales the observed kinematics and their errors by this set's
        projected masses, zeroes the first and last point of 1D orbit
        velocity histograms IN PLACE (mimicking ``triaxnnls_CRcut.f90``),
        transforms the orbit library into the observed parameterisation,
        and applies the CRcut if enabled (only has an effect for
        GaussHermite). Returns flat ``(obs_kins, obs_kins_err, orb_kins)``.

        The zeroing is undone before returning: the two edge slices are
        (n_orbs x n_apertures each - tiny), and handing the library back
        exactly as it was received lets callers reuse it instead of paying
        for a second full read, which is by far the most expensive operation
        in a model.
        """
        hist_dim = len(orb_veldist.y[0, ..., 0].shape)  # 1D or 2D vel hists
        obs_kins, obs_kins_err = tmp
        obs_kins = (obs_kins.T * prj_mass_i).T
        obs_kins_err = (obs_kins_err.T * prj_mass_i).T
        edges_restore = []
        if hist_dim == 1:  # Do we need this for proper motions (2d hists)?
            edges_restore.append((orb_veldist.y, orb_veldist.y[:, 0, :].copy(), orb_veldist.y[:, -1, :].copy()))
            orb_veldist.y[:, 0, :] = 0.0
            orb_veldist.y[:, -1, :] = 0.0
        # transform orblib to same parameterisation as observed kinematics
        orb_kins = kins.transform_orblib_to_observables(orb_veldist, self.settings)
        if self.CRcut:
            orb_kins = self.apply_CR_cut(kins, orb_veldist, orb_kins)
        # put the zeroed velocity-histogram edges back, so orblib is left
        # exactly as it was handed to us (see the note above)
        for hist_y, first, last in edges_restore:
            hist_y[:, 0, :] = first
            hist_y[:, -1, :] = last
        return np.ravel(obs_kins), np.ravel(obs_kins_err), orb_kins

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
        kins_data = kins.get_data()
        obs_mu_v = kins_data["v"]
        obs_sig_v = kins_data["sigma"]
        delta_v = np.abs(orb_mu_v - obs_mu_v)
        condition1 = np.abs(obs_mu_v) / obs_sig_v > 1.5
        condition2 = delta_v / obs_sig_v > 3.0
        condition3 = obs_mu_v * orb_mu_v < 0
        idx_cut = np.where(condition1 & condition2 & condition3)
        cut = np.zeros_like(orb_mu_v, dtype=bool)
        cut[idx_cut] = True
        naperture_cut = np.sum(cut, 1)
        # orbit 'j' is "bad" in naperture_cut[j] apertures
        # if an orbit is bad in 0 or 1 apertures, then we ignore this
        cut[naperture_cut < 1, :] = False
        # to cut an orbit, replace it's h1 by 3.0/dvhist(i)
        idx_cut = np.where(cut)
        dvhist = kins.hist_width / kins.hist_bins
        dvhist = np.max(dvhist)
        orb_gh[idx_cut[0], idx_cut[1], 0] = 3.0 / dvhist
        return orb_gh

    @staticmethod
    def kkt_violation(A, b, weights):
        r"""Optimality certificate for the NNLS solution.

        NNLS is convex, so the Karush-Kuhn-Tucker conditions are necessary
        **and sufficient**: a point satisfying them is a global optimum. With
        :math:`g = A^T(Aw-b)` they read, per orbit,

        - :math:`w_j > 0 \Rightarrow g_j = 0`
        - :math:`w_j = 0 \Rightarrow g_j \geq 0`   (adding weight cannot help)

        Two numbers are returned: the largest violation itself, and that
        violation scaled by ``||A_col|| * ||r||`` where ``r = Aw - b`` is the
        residual. By Cauchy-Schwarz :math:`|g_j| \leq \|A_{\cdot j}\|\,\|r\|`,
        so the scaled value lies in [0, 1] and measures how strongly the worst
        orbit is still aligned with the residual.

        The residual, not ``b``, is the correct denominator. Scaling by
        ``||b||`` was tried first and is unusable here: ``b[0] = 1e8`` makes the
        factor ~1e16, so a raw violation of 3.9e3 was reported as 3.9e-13 and a
        clearly unconverged solution looked optimal. That is the same failure
        mode as adelie's ``y_var``. ``||r||`` is a property of the solution and
        shrinks as it converges, which is what a convergence measure requires.

        Measured on solutions whose status is known independently:

        =========================== ========= ========== ==========
        solution                    chi2      raw        scaled
        =========================== ========= ========== ==========
        synthetic, exact            1.4e-32   1.1e-16    9.0e-09
        synthetic, scipy exact      9.2e-30   8.5e-15    2.8e-08
        synthetic, ALM unconverged  8.5e+00   4.0e+03    1.4e-05
        omega Cen, stored scipy     9.1e+09   7.1e+16    1.0e+00
        omega Cen, ALM              1.2e+04   4.5e+13    5.7e-03
        =========================== ========= ========== ==========

        The omega Cen scipy solution scores 1.0, the maximum: 27568 of its
        36000 orbits sit at zero while their gradients are negative, yet scipy
        reported convergence. chi2 cannot detect this, since comparing two
        solutions does not establish that either is optimal.

        Note what the table also shows: omega Cen's good solution (5.7e-03)
        scores WORSE than the synthetic unconverged one (1.4e-05). Being
        bounded in [0, 1] does not make the value comparable across problems,
        because what counts as small depends on the conditioning. Use it to
        rank solutions of one problem, and to catch gross failures near 1; do
        not read a fixed threshold as a convergence criterion.

        Parameters
        ----------
        A : array (n_constraints, n_orbits)
        b : array (n_constraints,)
        weights : array (n_orbits,)

        Returns
        -------
        tuple of float
            ``(scaled, raw)``. ``scaled`` is in [0, 1] and is comparable across
            problems; ``raw`` is dimensional and comparable only within one.
            Report both.

        References
        ----------
        Boyd, S. & Vandenberghe, L. 2004, Convex Optimization (Cambridge Univ.
            Press), Sect. 5.5.3 -- KKT conditions are sufficient, not only
            necessary, for a convex problem such as this one

        """
        resid = A @ weights - b
        grad = A.T @ resid
        viol = np.where(weights > 0, np.abs(grad), np.maximum(-grad, 0.0))
        raw = float(np.max(viol))
        # Cauchy-Schwarz denominator: |g_j| <= ||A_.j|| ||r||, so this is in
        # [0, 1]. An exactly-fitting solution has ||r|| -> 0 and the ratio is
        # then 0/0; guard it and report 0, which is the correct verdict.
        scale = np.linalg.norm(A, axis=0) * np.linalg.norm(resid)
        if not np.any(scale > 0):
            return 0.0, raw
        scale = np.where(scale > 0, scale, np.inf)
        return float(np.max(viol / scale)), raw

    @staticmethod
    def _build_augmented_X(A_rest, b_rest, sqrt_mu, dtype):
        """Build adelie's augmented, column-scaled matrix from A's body.

        Unit-L2 column scaling is an exact change of variable, since positive
        diagonal scaling preserves w >= 0, and the column norms span about 15
        orders of magnitude. The array is F-contiguous because coordinate
        descent accesses one column at a time.

        Bitwise identical to the inline construction this was extracted from
        (solve_adelie_alm, pre-2026-08-21); dev_tests/test_augmented_build.py
        pins that. Kept as a named step so the fused constructor
        (construct_adelie_matrix_and_rhs) shares the same finishing moves."""
        n_orbs = A_rest.shape[1]
        # Filled in place into an F-ordered buffer: np.vstack + `X / col_norm`
        # + np.asfortranarray would each allocate a full copy of X, so the
        # naive version peaks at 4x the matrix (~500 GiB for omega Cen).
        X = np.empty((A_rest.shape[0] + 1, n_orbs), dtype=dtype, order="F")
        X[0, :] = sqrt_mu
        X[1:, :] = A_rest
        col_norm, y = _scale_columns(X, b_rest, dtype)
        return X, col_norm, y

    @staticmethod
    def kkt_violation_augmented(row0_vec, b0, X_scaled, col_norm, resid_full, weights, mu):
        r"""Optimality certificate computed from the augmented matrix alone.

        Same semantics as :meth:`kkt_violation` — returns ``(scaled, raw)``
        with scaled in [0, 1] — but expressed through the fused augmented
        matrix so A never has to exist:

            A[1:, j]   = col_norm[j] * X_scaled[1:, j]
            grad[j]    = row0_vec[j]*r0 + col_norm[j]*(X_scaled[1:]^T r)[j]
            ||A_.j||^2 = row0_vec[j]^2 + col_norm[j]^2 - mu

        where ``r`` is the plain residual ``Aw - b`` aligned to A's rows
        (slot 0 = the total-mass row) and col_norm includes the penalty row.
        Two passes over X instead of three over A.

        Algebraically identical to kkt_violation(A, b, w); differs in rounding
        only. dev_tests/test_surrogate_chi2_kkt.py checks agreement to rtol
        1e-10 and the exact-fit / degenerate-column guards.
        """
        r0 = float(row0_vec @ weights) - float(b0)
        r_rest = np.asarray(resid_full, dtype=np.float64).ravel()[1:]
        resid = np.concatenate(([r0], r_rest))
        grad = row0_vec.astype(np.float64) * r0 + col_norm.astype(np.float64) * (X_scaled[1:, :].T @ r_rest)
        viol = np.where(weights > 0, np.abs(grad), np.maximum(-grad, 0.0))
        raw = float(np.max(viol))
        # Cauchy-Schwarz denominator as in kkt_violation, rebuilt through the
        # identities above; round-off can make the expression dip below zero,
        # hence the clip. An exactly-fitting solution reports 0 like the stock
        # version's 0/0 guard.
        col_sq = np.maximum(
            col_norm.astype(np.float64) ** 2 - mu + row0_vec.astype(np.float64) ** 2,
            0.0,
        )
        scale = np.sqrt(col_sq) * np.linalg.norm(resid)
        if not np.any(scale > 0):
            return 0.0, raw
        scale = np.where(scale > 0, scale, np.inf)
        return float(np.max(viol / scale)), raw

    def solve_adelie_alm(self, problem):
        r"""Solve the NNLS problem with adelie BVLS + an augmented Lagrangian.

        Consumes an :class:`AdelieProblem` built by
        :meth:`construct_adelie_matrix_and_rhs` (or
        :meth:`_build_augmented_X` from a classic A, b pair) and returns
        ``(best_w, resid_full)``: the best weights and the plain residual
        ``A @ w - b`` aligned to A's rows (slot 0 = total-mass row), so the
        caller can compute chi2 without A.

        **Why an augmented Lagrangian.** Row 0 enforces
        :math:`\sum_j w_j = M_{tot}` with ``econ[0] = 1e-8``, so it contributes
        :math:`\tfrac{1}{2}\,10^{16}(\mathbf{1}^Tw-1)^2`. That is the
        :math:`\mu\to\infty` limit of a quadratic penalty, i.e. a hard equality
        constraint imposed by brute force. On omega Cen the row raises the
        condition number to 5e22. It also inflates adelie's convergence
        threshold: the stopping test is ``convg_measure <= tol * y_var`` with
        ``y_var = ||b||^2/n``, and ``b[0] = 1e8`` gives ``y_var = 1.04e12``, so
        the threshold becomes ~1e5 instead of ~1e-7. adelie then reports
        convergence after 2 iterations at ``sum(w) = 0``.

        ALM imposes the same constraint with a moderate ``mu`` by carrying a
        multiplier, so the matrix passed to adelie stays well scaled::

            row     = sqrt(mu) * 1',   target = sqrt(mu) * (1 + lam/mu)
            lam    <- lam - mu * (1'w - 1)

        A fixed penalty alone would leave a biased optimum. Updating the target
        through ``lam`` removes that bias without raising ``mu``: on omega Cen
        the constraint gap reaches 6e-11 with ``mu = 1e7``.

        Three implementation details are needed for this to be practical:

        1. ``mu`` is fixed. Raising it when convergence stalls restores the
           original scaling problem; one such attempt reached ``mu = 5.3e11``
           and stalled at a gap of 2e-4.
        2. With ``mu`` fixed the augmented matrix does not change between
           iterations (only the scalar target does), so it is built and
           column-scaled once.
        3. Each subproblem is warm-started from the previous adelie state.
           Consecutive subproblems differ by one number and converge in 1 to 2
           inner iterations.

        4. chi2 for the iterate is taken from ``state.resid`` rather than
           recomputed as ``A @ w - b``. The column scaling is exact, so
           ``X[1:] @ beta == A[1:] @ w`` and adelie's residual already holds
           everything except row 0 of A, which is one dot product over the
           orbits. Evaluating ``A @ w`` instead costs a full pass over the
           matrix (125 GiB for omega Cen) on every multiplier update.

        Items 2 and 3 reduced the per-iteration cost from 18.2 s to 0.6 s on
        NGC6278, which is what makes several hundred multiplier updates
        affordable. Item 4 matters on omega Cen rather than NGC6278: it is
        proportional to the size of the matrix, not to the solve.

        The inner solves are inexact, so the multiplier oscillates about its
        fixed point and the final iterate is arbitrary within that spread. The
        iterate with the lowest chi2 is returned instead.

        .. warning::
           **Known limitation: ``adelie_mu`` is an absolute constant but should
           be relative to the matrix scale.** The default 1e7 was validated on
           omega Cen and NGC6278, where the non-total-mass rows have column
           norms of order 1e4 to 1e5, so ``sqrt(mu) ~ 3e3`` is comparable to
           them. On a synthetic problem whose other rows are of order 1, the
           same mu leaves chi2 = 8.5 where 0 is attainable, while mu = 1e2
           reaches 1.6e-4. A galaxy whose constraint magnitudes differ from
           omega Cen's may therefore sit off the validated plateau, and the
           failure is silent. Until mu is set from the data, check the logged
           chi2 against a scipy solve on any new dataset.

        Parameters
        ----------
        problem : AdelieProblem

        Returns
        -------
        tuple (array, array)
            ``(best_w, resid_full)``: orbit weights and the plain residual
            aligned to A's rows.

        References
        ----------
        Hestenes, M. R. 1969, J. Optim. Theory Appl., 4, 303 -- augmented
            Lagrangian method (independently of Powell)
        Powell, M. J. D. 1969, in Optimization, ed. R. Fletcher (Academic
            Press), 283
        Nocedal, J. & Wright, S. J. 2006, Numerical Optimization, 2nd edn.
            (Springer), Chap. 17 -- standard treatment, including why the
            multiplier update removes the bias of a fixed finite penalty
        Yang, J. & Hastie, T. 2024a, arXiv:2410.03014 -- the BVLS coordinate
            descent solver called here
        Yang, J. & Hastie, T. 2024b, arXiv:2405.08631 -- the adelie package

        """
        n_threads = int(os.environ.get("OMP_NUM_THREADS", os.cpu_count() or 1))
        mu = self.adelie_mu
        sqrt_mu = np.sqrt(mu)
        X, col_norm, y = problem.X, problem.col_norm, problem.y
        dtype = X.dtype
        n_orbs = X.shape[1]

        lower = np.zeros(n_orbs, dtype=dtype)
        upper = np.full(n_orbs, np.inf, dtype=dtype)
        # adelie's bvls() infers a state dtype from X but defaults `weights`
        # to np.full(n, 1/n) (always float64); when X is float32 that
        # mismatch makes an internal `np.array(..., copy=False)` raise under
        # numpy>=2.0. Passing weights explicitly in X's dtype avoids it.
        weights_arr = np.full(X.shape[0], 1 / X.shape[0], dtype=dtype)

        lam = 0.0
        state = None
        best_chi2, best_w, best_it = np.inf, None, -1
        for it in range(self.adelie_alm_iters):
            y[0] = sqrt_mu * (1.0 + lam / mu)
            state = _adelie_solver.bvls(
                X,
                np.ascontiguousarray(y),
                lower,
                upper,
                weights=weights_arr,
                n_threads=n_threads,
                tol=self.adelie_tol,
                max_iters=int(2e5),
                warm_start=state,
            )
            w = (np.asarray(state.beta).ravel() / col_norm).astype(np.float64)
            gap = float(w.sum() - 1.0)
            lam -= mu * gap
            # chi2 without touching A at all. X[1:] is A[1:] with unit-L2
            # column scaling and w = beta/col_norm, so X[1:] @ beta is exactly
            # A[1:] @ w, and adelie already returns resid = y - X @ beta from
            # the solve we just did. Only row 0 of A needs evaluating here,
            # because X replaces it with the ALM penalty row. This drops a full
            # pass over A (~125 GiB for omega Cen) from every iteration.
            # The astype is a no-op at the default nnls_dtype and a deliberate
            # 3 MB upcast at float32, where summing 371212 float32 terms would
            # lose precision that the old float64 accumulation had.
            resid = np.asarray(state.resid).ravel()[1:]
            resid = resid.astype(np.float64, copy=False)
            # row0_vec is bitwise A[0], so this equals the old float(A[0]@w)
            row0 = float(problem.row0_vec @ w) - problem.b0
            chi2 = row0 * row0 + float(resid @ resid)
            if chi2 < best_chi2:
                best_chi2, best_w, best_it = chi2, w.copy(), it
            if abs(gap) < self.adelie_gap_tol:
                break

        self.logger.info(
            f"adelie ALM: {it + 1} iterations, mu={mu:.1e}, "
            f"final gap={gap:.2e}, best iterate {best_it}, "
            f"chi2={best_chi2:.4f}, sum(w)={best_w.sum():.10f}"
        )
        # plain residual at the returned weights, aligned to A's rows, built
        # from X: rows 1.. via the exact column scaling (round-off-level
        # difference vs a gemv over A), row 0 via row0_vec. Serves both the
        # surrogate KKT below and solve()'s final chi2_vector.
        r0 = float(problem.row0_vec @ best_w) - problem.b0
        r_rest = y[1:].astype(np.float64) - X[1:, :] @ (col_norm.astype(np.float64) * best_w)
        resid_full = np.concatenate(([r0], r_rest))
        kkt, kkt_raw = self.kkt_violation_augmented(problem.row0_vec, problem.b0, X, col_norm, resid_full, best_w, mu)
        self.logger.info(f"adelie ALM: KKT violation scaled={kkt:.3e} (in [0,1]), raw={kkt_raw:.3e}")
        # A fixed threshold cannot separate converged from unconverged across
        # datasets: omega Cen's good solution scores 5.7e-03 while a known
        # UNCONVERGED synthetic one scores 1.4e-05. What counts as small is
        # problem dependent. 0.1 therefore only catches unambiguous failures
        # (the stored omega Cen scipy solution scores 1.0, the maximum).
        # For a real convergence check, compare chi2 against a scipy solve.
        if kkt > 0.1:
            self.logger.warning(
                f"adelie ALM: scaled KKT violation {kkt:.3e} is close to the "
                "maximum of 1 - the solution is far from optimal. Check "
                "adelie_mu against a scipy solve on this dataset."
            )
        return best_w, resid_full

    def solve(self, orblib, ignore_existing_weights=False):
        """Solve for orbit weights

        **Note:** the returned chi2 values are not the same as
        ``LegacyWeightSolver.read_chi2`` - see the docstring for more info

        Apart from weight solving, the attributes ``orblib.intrinsic_masses``,
        ``orblib.projected_masses``, and ``orblib.vel_histograms`` are set
        via calling ``orblib.read_vel_histograms()``.

        Parameters
        ----------
        orblib : dyn.OrbitLibrary
        ignore_existing_weights : bool, optional
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
        self.logger.info(f"Using WeightSolver: {__class__.__name__}/{self.nnls_solver}")
        if (not ignore_existing_weights) and self.weight_file_exists():
            results = ascii.read(self.weight_file, format="ecsv")
            self.logger.info("NNLS solution read from existing output")
            weights = results["weights"]
            chi2_tot = results.meta["chi2_tot"]
            chi2_kin = results.meta["chi2_kin"]
            chi2_kinmap = results.meta["chi2_kinmap"]
        else:
            # On the fused+streaming path the constructor drives per-set reads
            # itself (freeing each set after use), so no eager full read here.
            adelie_streaming = self.nnls_solver == "adelie" and self.stream_reads
            # chi2_kinmap needs the library exactly as assembly received it.
            # Assembly restores the edge bins it zeroes, but streaming frees
            # each set's histograms and float32 rewrites hist.y in place.
            orblib_reusable = not adelie_streaming and self.nnls_dtype != np.float32
            if not adelie_streaming:
                orblib.read_vel_histograms()  # sets orblib.vel_histograms,
                # orblib.intrinsic_masses, and
                # orblib.projected_masses
                # downcast the retained orbit-library data before building the
                # NNLS matrix (no-op unless nnls_dtype is float32)
                _downcast_orblib(orblib, self.nnls_dtype)
            if self.nnls_solver == "adelie":
                # The fused constructor assembles adelie's augmented matrix X
                # directly from the orbit library - A is never built on this
                # path, which halves the resident set during the solve.
                if not _ADELIE_AVAILABLE:
                    text = "nnls_solver 'adelie' is not installed. Run: pip install adelie"
                    self.logger.error(text)
                    raise ImportError(text)
                try:
                    problem = self.construct_adelie_matrix_and_rhs(orblib)
                    weights, resid_full = self.solve_adelie_alm(problem)
                except Exception as e:
                    txt = (
                        f"Orblib {orblib.mod_dir}, ml={orblib.parset['ml']}"
                        f": adelie ALM solver error occured: {e} All weights "
                        "and chi2 set to nan. Consider trying scipy."
                    )
                    self.logger.warning(txt)
                    weights = np.full(orblib.n_orbs, np.nan)
            elif self.nnls_solver == "scipy":
                # scipy/cvxopt keep the classic A-based path. The adelie
                # branch above builds its augmented matrix directly instead.
                A, b = self.construct_nnls_matrix_and_rhs(orblib)
                # Normalize the data: building it on the adelie path would
                # hold a second full copy of A (~125 GiB for omega Cen).
                A_max = np.max(np.abs(A), axis=0)
                A_normalized = A / A_max
                b_max = np.max(np.abs(b))
                b_normalized = b / b_max
                try:
                    # Solve the NNLS problem with normalized data
                    x_normalized, rnorm = optimize.nnls(A_normalized, b_normalized)
                    # Scale the solution back to the original scale
                    weights = x_normalized * b_max / A_max

                except Exception as e:
                    txt = (
                        f"Orblib {orblib.mod_dir}, ml={orblib.parset['ml']}"
                        f": SciPy solver error occured: {e} All weights "
                        "and chi2 set to nan. Consider trying cvxopt."
                    )
                    self.logger.warning(txt)
                    weights = np.full(A.shape[1], np.nan)
            elif self.nnls_solver == "cvxopt":
                A, b = self.construct_nnls_matrix_and_rhs(orblib)
                A_max = np.max(np.abs(A), axis=0)
                A_normalized = A / A_max
                b_max = np.max(np.abs(b))
                b_normalized = b / b_max
                try:
                    P = np.dot(A_normalized.T, A_normalized)
                    q = -1.0 * np.dot(A_normalized.T, b_normalized)
                    solver = CvxoptNonNegSolver(P, q)
                    x_normalized = solver.beta
                    # Scale the solution back to the original scale
                    weights = x_normalized * b_max / A_max
                except Exception as e:
                    txt = (
                        f"Orblib {orblib.mod_dir}, ml={orblib.parset['ml']}"
                        f": CVXOPT solver error occured: {e} All weights "
                        "and chi2 set to nan. Consider trying scipy."
                    )
                    self.logger.warning(txt)
                    weights = np.full(A.shape[1], np.nan)
            else:
                text = "Unknown nnls_solver"
                self.logger.error(text)
                raise ValueError(text)
            if not np.isnan(weights[0]):
                # calculate chi2s
                if self.nnls_solver == "adelie":
                    # chi2 without A: the plain residual at the returned
                    # weights, aligned to A's rows, plus the total-mass row
                    # term. Algebraically identical to (A@w - b)**2; differs
                    # only in gemv rounding (quantified in the perf notes).
                    row0_resid = float(problem.row0_vec @ weights) - problem.b0
                    chi2_vector = chi2_vector_from_residuals(resid_full, row0_resid * row0_resid)
                    # X is dead from here on, and chi2_kinmap below may read a
                    # fresh orbit library (~165 GiB at production scale). Drop
                    # X first so the two peaks do not stack - glibc unmaps
                    # numpy's buffers at `del`, measured at 100% reclaim.
                    del problem, resid_full
                else:
                    chi2_vector = (np.dot(A, weights) - b) ** 2.0
                chi2_tot = np.sum(chi2_vector)
                chi2_kin = np.sum(chi2_vector[self.n_mass_constraints :])
                chi2_kinmap = self.chi2_kinmap(weights, orblib=orblib if orblib_reusable else None)
                # save the output
                results = table.Table()
                results["weights"] = weights
                # add chi2 to meta data
                results.meta = {
                    "chi2_tot": chi2_tot,
                    "chi2_kin": chi2_kin,
                    "chi2_kinmap": chi2_kinmap,
                }
                results.write(self.weight_file, format="ascii.ecsv", overwrite=True)
                self.logger.info("NNLS problem solved and chi2 calculated.")
            else:
                chi2_tot = chi2_kin = chi2_kinmap = np.nan
            # delete existing .yaml files and copy current config file
            # into model directory
            self.config.copy_config_file(self.direc_with_ml)
        return weights, chi2_tot, chi2_kin, chi2_kinmap


class CvxoptNonNegSolver:
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
        self.success = sol["status"] == "optimal"
        self.beta = np.squeeze(np.array(sol["x"]))


# end
