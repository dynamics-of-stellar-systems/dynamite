import os
import copy
import glob
import difflib
import logging
import shutil
import numpy as np
from astropy import table
from astropy.io import ascii

from dynamite import weight_solvers as ws
from dynamite import orblib as dyn_orblib

class AllModels(object):
    """All models which have been run so far

    The main attribute ``self.table`` is an Astropy table holding all models run
    so far.

    Parameters
    ----------
    config : a ``dyn.config_reader.Configuration`` object
    from_file : bool
        whether to create this object from a saved `all_models.ecsv` file

    """
    def __init__(self, config=None, from_file=True):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        if config is None:
            text = f'{__class__.__name__} needs configuration object, ' \
                   'None provided.'
            self.logger.error(text)
            raise ValueError(text)
        self.config = config
        self.system = config.system
        self.set_filename(config.settings.io_settings['all_models_file'])
        self.make_empty_table()
        if from_file and os.path.isfile(self.filename):
            self.logger.info('Previous models have been found: '
                             f'Reading {self.filename} into '
                             f'{__class__.__name__}.table')
            self.read_model_table()
        else:
            self.logger.info(f'No previous models (file {self.filename}) '
                             'have been found: Made '
                             f'an empty table in {__class__.__name__}.table')

    def set_filename(self, filename):
        """Set the name (including path) for this model

        Parameters
        ----------
        filename : string
            name for this model

        Returns
        -------
        None
            sets ``self.filename``

        """
        outdir = self.config.settings.io_settings['output_directory']
        filename = f'{outdir}{filename}'
        self.filename = filename

    def make_empty_table(self):
        """Make an empty Astropy table

        Returns
        -------
        None
            sets ``self.table``

        """
        names = self.config.parspace.par_names.copy()
        dtype = [float for n in names]
        # add the columns from legacy version
        names += ['chi2', 'kinchi2', 'kinmapchi2', 'time_modified']
        dtype += [float, float, float, 'U256']
        # add extra columns
        names += ['orblib_done', 'weights_done', 'all_done']
        dtype += [bool, bool, bool]
        # which_iter will record which iteration of parameters a model came from
        names.append('which_iter')
        dtype.append(int)
        # directory will be the model directory name in the models/ directory
        names.append('directory')
        dtype.append('U256')
        self.table = table.Table(names=names, dtype=dtype)

    def read_model_table(self):
        """Read table from file ``self.filename``

        Returns
        -------
        None
            sets ``self.table``

        """
        table_read = ascii.read(self.filename)
        self.table = table.vstack((self.table, table_read),
                                  join_type='outer',
                                  metadata_conflicts='error')
        self.logger.debug(f'{len(self.table)} models read '
                          f'from file {self.filename}')

    def update_model_table(self):
        """all_models table update: fix incomplete models, add kinmapchi2.

        Dealing with incomplete models:
        Models with all_done==False but an existing model_done_staging.ecsv
        will be updated in the table and the staging file will be deleted.
        Models with all_done==False and no existing orblib will be deleted
        from the table and their model directory will be deleted, too.
        The configuration setting reattempt_failures determines how partially
        completed models with all_done==False but existing orblibs are treated:
        If reattempt_failures==True, their orblib_done will be set to True
        and later the ModelIterator will execute the weight solving
        based on the existing orblibs.
        If reattempt_failures==False, the model and its directory will be
        deleted.
        Note that orbit libraries on disk will not be deleted as they
        may be in use by other models.

        Up to DYNAMITE 3.0 there was no kinmapchi2 column in the all_models
        table. If possible (data exists on disk), calculate and add the values,
        otherwise set to np.nan.

        Returns
        -------
        None
            sets ``self.table``

        """
        table_modified = False
        for i, row in enumerate(self.table):
            if not row['all_done']:
                table_modified = True
                mod = self.get_model_from_row(i)
                staging_filename = mod.directory+'model_done_staging.ecsv'
                f_root = mod.directory_noml + 'datfil/'
                check = os.path.isfile(f_root + 'orblib.dat.bz2') \
                        and os.path.isfile(f_root + 'orblibbox.dat.bz2')
                if not check:
                    check = os.path.isfile(f_root + 'orblib_qgrid.dat.bz2') \
                     and os.path.isfile(f_root + 'orblib_losvd_hist.dat.bz2') \
                     and os.path.isfile(f_root + 'orblibbox_qgrid.dat.bz2') \
                     and os.path.isfile(f_root + 'orblibbox_losvd_hist.dat.bz2')
                if os.path.isfile(staging_filename):
                    # the model has completed but was not entered in the table
                    staging_file = ascii.read(staging_filename)
                    self.table[i] = staging_file[0]
                    self.logger.info(f'Staging file {staging_filename} '
                                f'used to update {__class__.__name__}.table.')
                    os.remove(staging_filename)
                    self.logger.debug(
                        f'Staging file {staging_filename} deleted.')
                elif check:
                    self.logger.debug(f'Row {i}: orblibs were computed '
                                      'but not weights.')
                    self.table[i]['orblib_done'] = True
                else:
                    self.logger.debug(f'Row {i}: neither orblibs nor '
                                      'weights were completed.')
        # collect failed models to delete (both their directory and table entry)
        to_delete = []
        # if we will reattempt weight solving, only delete models with no orblib
        if self.config.settings.weight_solver_settings['reattempt_failures']:
            for i, row in enumerate(self.table):
                if (not row['orblib_done']) and (not row['all_done']):
                    to_delete.append(i)
                    self.logger.info('No orblibs calculated for model in '
                                     f'{row["directory"]} - removing row {i}.')
        # otherwise delete any model which is not `all_done`
        else:
            for i, row in enumerate(self.table):
                if not row['all_done']:
                    to_delete.append(i)
                    self.logger.info('No finished model found in '
                                     f'{row["directory"]} - removing row {i}.')
        # do the deletion
        # note: only models without orblibs are deleted, so we delete the
        # entire orblibs' directories
        if len(to_delete)>0:
            cwd = os.getcwd()
            os.chdir(self.config.settings.io_settings['model_directory'])
            dirs_to_delete = set(
                                 [d[:d[:-1].rindex('/')+1]
                                 for d in self.table[to_delete]['directory']]
                                )
            self.logger.info(f'Will try to remove {len(dirs_to_delete)} '
                             'unique orblibs.')
            for directory in dirs_to_delete:
                try:
                    # Only remove orblib directories that are not used by
                    # already completed models
                    orblibs_keep = [d[:d[:-1].rindex('/')+1] for d in
                     self.table[np.where(self.table['all_done'])]['directory']]
                    if directory in set(orblibs_keep):
                        self.logger.info(f'Orblib directory {directory} in '
                                         'use by existing model - untouched.')
                    else:
                        shutil.rmtree(directory)
                        self.logger.info(f'Orblib directory {directory} '
                                         'removed.')
                except:
                    self.logger.warning(f'Cannot remove orblib in {directory},'
                        ' perhaps it has already been removed before.')
            os.chdir(cwd)
            self.table.remove_rows(to_delete)
        # Up to DYNAMITE 3.0 there was no kinmapchi2 column -> retrofit.
        if isinstance(self.table['kinmapchi2'], table.column.MaskedColumn):
            table_modified = True
            self.retrofit_kinmapchi2()
        # If the table has been modified, save it.
        if table_modified:
            self.save()
            self.logger.info('all_models table updated and saved.')
        else:
            self.logger.info('No all_models table update required.')

    def retrofit_kinmapchi2(self):
        """Calculates kinmapchi2 for DYNAMITE legacy tables if possible.

        Returns
        -------
        None.
            updates ``self.table``
        """
        which_chi2 = 'kinmapchi2'
        self.logger.info('Legacy all_models table read, updating '
                         f'{which_chi2} column...')
        # self.table[which_chi2] = np.nan
        for row_id, row in enumerate(self.table):
            if row['orblib_done'] and row['weights_done']:
                # both orblib_done==True and weights_done==True indicates
                # that data for kinmapchi2 is on the disk -> calculate
                # kinmapchi2
                mod = self.get_model_from_row(row_id)
                ws_type = self.config.settings.weight_solver_settings['type']
                weight_solver = getattr(ws, ws_type)(
                                        config=self.config,
                                        directory_with_ml=mod.directory)
                orblib = mod.get_orblib()
                _, _, _, row[which_chi2] = weight_solver.solve(orblib)
                self.logger.info(f'Model {row_id}: {which_chi2} = '
                                 f'{row[which_chi2]}')
            else:
                row[which_chi2] = np.nan
                self.logger.warning(f'Model {row_id}: cannot update '
                                    f'{which_chi2} - data deleted?')

    def read_legacy_chi2_file(self, legacy_filename):
        """
        Read the `legacy` AKA schwpy format of chi2 files

        Taken from schw_basics.py, reads in legacy files named similar to
        griddata/_chi2.cat

        Parameters
        ----------
        legacy_filename: string
            the legacy_filename (probably griddata/_chi2.cat)
        """
        ### read the header
        head1 = np.genfromtxt(legacy_filename, max_rows=1)
        Nf = int(head1[0])
        npar = int(head1[1])
        ### read the main matrix
        mpar = np.genfromtxt(legacy_filename,
                             max_rows=Nf,
                             skip_header=1)
        mtest = np.genfromtxt(legacy_filename,
                              max_rows=1,
                              skip_header=Nf + 1)
        len_mtest = len(np.atleast_1d(mtest))
        ### read the last modification time
        if np.mod(Nf, len_mtest) > 0:
            mlast = 1  # add 1 to line counts for the next step, reading the fls
            mtime1 = np.genfromtxt(legacy_filename,
                                   max_rows=int(Nf / len_mtest),
                                   skip_header=Nf + 1)
            mtime2 = np.genfromtxt(legacy_filename,
                                   max_rows=1,
                                   skip_header=Nf + int(Nf / len(mtest)) + 1)
            mtime = np.hstack((np.ravel(mtime1), np.ravel(mtime2)))
        else:
            mlast = 0
            mtime = np.ravel(np.genfromtxt(legacy_filename,
                                           max_rows=int(Nf / len_mtest),
                                           skip_header=Nf + 1))
        ### read the file paths
        fls = np.genfromtxt(legacy_filename,
                            max_rows=Nf,
                            skip_header=Nf + int(Nf / len_mtest) + 1 + mlast,
                            dtype=str)
        self.logger.debug(f'Legacy chi2 file {legacy_filename} read')
        return Nf, npar, mpar.T, mtime.T, fls.T

    def convert_legacy_chi2_file(self, legacy_filename=None):
        """
        Convert `legacy` format of chi2 files

        `legacy` AKA schwpy format were likely called ```griddata/_chi2.cat``.

        Parameters
        ----------
        legacy_filename: string
            the legacy_filename (probably griddata/_chi2.cat)
        """
        # TODO: (maybe...?)
        # make more general if legacy parameters have different names
        Nf, npar, mpar, mtime, fls = self.read_legacy_chi2_file(legacy_filename)
        # legacy parameter names are fixed: bh, dc, etc...
        mods = table.Table()
        mods['bh'] = mpar[0,:]
        mods['dc'] = mpar[1,:]
        mods['f'] = mpar[2,:]
        mods['q'] = mpar[3,:]
        mods['p'] = mpar[4,:]
        mods['u'] = mpar[5,:]
        mods['ml'] = mpar[6,:]
        mods['chi2'] = mpar[7,:]
        mods['kinchi2'] = mpar[8,:]
        mods['time_modified'] = mtime
        # currently fls are strings "{param_directory}/{ml_directory}/nn"
        # cleaner to have just the directry name instead...?
        direcs = [file0[:-2] for file0 in fls]
        mods['directory'] = table.Column(direcs, dtype='<U100')
        # add any extra columns we want to have in the table
        mods['ics_done'] = True
        mods['orblib_done'] = True
        mods['weights_done'] = True
        mods['all_done'] = True
        # which_iter will record which iteration of parameters a model came from
        mods['which_iter'] = 0  # was not recorded for schwpy so set to 0
        mods.write(self.filename, format='ascii.ecsv')
        self.logger.debug(f'Legacy chi2 file {legacy_filename} converted ' + \
                          f'to {self.filename}')

    def get_parset_from_row(self, row_id):
        """Get a parset set given a table row

        Parameters
        ----------
        row_id : int
            which row

        Returns
        -------
        list
            a list of ``dyn.parspace.Parameter`` objects

        """
        parset = self.table[row_id][self.config.parspace.par_names]
        return parset

    def get_model_from_parset(self, parset):
        """Get the ``Model`` from a parset

        Parameters
        ----------
        parset : list
            a list of ``dyn.parspace.Parameter`` objects

        Returns
        -------
        a ``dyn.model.Model`` object

        """
        for idx, row in enumerate(self.table[self.config.parspace.par_names]):
            if np.allclose(tuple(parset), tuple(row)):
                mod = self.get_model_from_row(idx)
                break
        else:
            text = f'parset not in all_models table. parset={parset}, ' \
                   f'all_models table: {self.table}'
            self.logger.error(text)
            raise ValueError(text)
        return mod

    def get_model_from_directory(self, directory):
        """Get the ``Model`` from a model directory

        Parameters
        ----------
        directory : str
            The directory string needs to start with the output directory
            defined in ``config.settings.io_settings['output_directory']``

        Raises
        ------
        ValueError
            If the directory does not exist in the all_models table.

        Returns
        -------
        mod : a ``dyn.model.Model`` object

        """
        for idx, dir_table in enumerate(self.table['directory']):
            if self.config.settings.io_settings['output_directory'] \
                                        + 'models/' + dir_table == directory:
                mod = self.get_model_from_row(idx)
                break
        else:
            text = f'Directory {directory} not in all_models table. ' \
                   f'all_models table: {self.table}'
            self.logger.error(text)
            raise ValueError(text)
        return mod

    def get_model_from_row(self, row_id):
        """Get a ``Model`` given a table row

        Parameters
        ----------
        row_id : int
            which row

        Returns
        -------
        a ``dyn.model.Model`` object

        """
        parset = self.get_parset_from_row(row_id)
        mod = Model(config=self.config,
                    parset=parset,
                    directory=self.table['directory'][row_id])
        return mod

    def get_row_from_model(self, model=None):
        """Get the table row for a ``Model``

        Parameters
        ----------
        model : a ``Model`` object

        Raises
        ------
        ValueError
            If the model is not found in ``self.table``.

        Returns
        -------
        row_id : int
            The ``self.table`` row index of model.

        """
        row_comp = tuple(model.parset[self.config.parspace.par_names])
        for row_id, row in enumerate(
                                self.table[self.config.parspace.par_names]):
            if np.allclose(row_comp, tuple(row)):
                break
        else:
            text = 'Cannot find model with parset ' \
                   f'{row_comp} in all_models table.'
            self.logger.error(text)
            raise ValueError(text)
        return row_id

    def get_ml_of_original_orblib(self, model_id):
        """Get ``ml`` of model number model_id's original orblib

        The original ``ml`` is required to rescale orbit libraries for rescaled
        potentials. This method searches ``self.table``, the all_models table.

        Parameters
        ----------
        model_id : int
            The ``self.table`` row index of the model.

        Raises
        ------
        ValueError
            If the ``ml`` parameter is not in the parameter space or the
            model's orblib cannot be found.

        Returns
        -------
        ml_orblib : float
            the original ``ml``

        """
        orblib_parameters = self.config.parspace.par_names[:]
        ml = 'ml'
        try:
            orblib_parameters.remove(ml)
        except:
            self.logger.error(f"Parameter '{ml}' not found - check "
                              "implementation")
            raise
        row_comp = tuple(self.table[orblib_parameters][model_id])
        model_dir = self.table['directory'][model_id]
        for row_id, row in enumerate( \
                                 self.table[orblib_parameters][:model_id+1]):
            if np.allclose(row_comp, tuple(row)):
                ml_orblib = self.table['ml'][row_id]
                ml_orblib_dir = self.table['directory'][row_id]
                self.logger.debug(f'Orblib of model {model_dir} has original '
                                  f'ml value of {ml_orblib} '
                                  f'(model {ml_orblib_dir}).')
                break
        else:
            text = f'Cannot find orblib of model {model_dir} in ' \
                   'all_models table.'
            self.logger.error(text)
            raise ValueError(text)
        return ml_orblib

    def get_model_velocity_scaling_factor(self, model_id=None, model=None):
        """Get the model's velocity scaling factor

        Returns sqrt(model_ml/original_orblib_ml).
        The model can be either given by its row id in ``self.table`` or
        as a ``Model`` object. Note that the parameters model_id and model
        are mutually exclusive.

        Parameters
        ----------
        model_id : int compatible
            The model's row id in ``self.table``.
        model : a ``Model`` object

        Raises
        ------
        ValueError
            If not exactly one of model_id and model are supplied.

        Returns
        -------
        scaling_factor : float
            The model's velocity scaling factor
            sqrt(model_ml/original_orblib_ml).

        """
        if model_id is None and isinstance(model, Model):
            model_id = self.get_row_from_model(model)
        elif not (model_id == int(model_id) and model is None):
            text = 'Need to pass either model_id (int) or model (Model).'
            self.logger.error(text)
            raise ValueError(text)
        ml_orblib = self.get_ml_of_original_orblib(model_id)
        scaling_factor = np.sqrt(self.table['ml'][model_id]/ml_orblib)
        return scaling_factor

    def save(self):
        """Save the all_models table

        """
        self.table.write(self.filename, format='ascii.ecsv', overwrite=True)
        self.logger.debug(f'Model table written to file {self.filename}')

    def get_best_n_models(self, n=10, which_chi2=None):
        """Get the best n models or all but the n best models so far

        Parameters
        ----------
        n : int, optional
            How many models to get. If negative, all models except the
            n best models will be returned. The default is 10.
        which_chi2 : str, optional
            Which chi2 is used for determining the best models. If None, the
            setting from the configuration file will be used.
            The default is None.

        Raises
        ------
        ValueError
            If which_chi2 is neither None nor a valid chi2 type.

        Returns
        -------
        a new ``astropy.table`` object holding the best n models, sorted by
        which_chi2

        """
        which_chi2 = self.config.validate_chi2(which_chi2)
        table = copy.deepcopy(self.table)
        table.sort(which_chi2)
        if n>=0:
            table = table[:n]
        else:
            table = table[-n:]
        return table

    def get_best_n_models_idx(self, n=10, which_chi2=None):
        """Get the indices of the best n models so far

        Parameters
        ----------
        which_chi2 : str, optional
            Which chi2 is used for determining the best models. If None, the
            setting from the configuration file will be used.
            The default is None.

        Returns
        -------
        list of int
            indices in the all_models table of the n best model so far, sorted
            by which_chi2

        """
        which_chi2 = self.config.validate_chi2(which_chi2)
        return list(self.table.argsort(keys=which_chi2)[:n])

    def get_mods_within_chi2_thresh(self, which_chi2=None, delta=None):
        """Get models within or outside a delta threshold of the best

        Parameters
        ----------
        which_chi2 : str, optional
            Which chi2 is used for determining the best models. If None, the
            setting from the configuration file will be used.
            The default is None.
        delta : float, optional
            The threshold value. Models with chi2 values differing
            from the opimum by at most delta will be returned. If None,
            models within 10% of the optimal value will be returned. If
            delta is negative, models that are NOT within a delta
            threshold of the best are returned.
            The default is to return models within 10% of the best.

        Raises
        ------
        ValueError
            If which_chi2 is neither None nor a valid chi2 type.

        Returns
        -------
        a new ``astropy.table`` object holding the ''delta-best'' models
        (if delta >= 0) or holding all but the ''delta-best'' models
        (if delta < 0), respectively.

        """
        which_chi2 = self.config.validate_chi2(which_chi2)
        chi2_min = np.nanmin(self.table[which_chi2])
        if delta is None:
            delta = chi2_min * 0.1
        if delta >= 0:
            models = self.table[self.table[which_chi2] <= chi2_min+delta]
        else:
            models = self.table[self.table[which_chi2] > chi2_min-delta]
        return models

    def make_best_models_table(self,
                               which_chi2=None,
                               n=None,
                               delta=None,
                               filename=None):
        """Make a table of the best models and save it to disk

        Parameters
        ----------
        which_chi2 : str, optional
            Which chi2 is used for determining the best models. If None, the
            setting from the configuration file will be used.
            The default is None.
        n : int, optional
            How many models to get. If None, n will be ignored.
            Default: if delta is specified, the default is none; if delta
            is None, the default is 10.
        delta : float, optional
            The threshold value. Models with chi2 values differing
            from the opimum by at most delta will be returned. If None,
            delta will be ignored. The default is None.
        filename : str, optional
            File name of the best models table. The file is written into the
            output directory specified in the config file. If None, the name
            is the same as for the all models table but with '_best' added
            to the base file name. If the file already exists, a warning
            is logged, the existing file is backed up ('_backup' added),
            and then overwritten. The default is None.

        Raises
        ------
        ValueError
            If both n and delta are specified (i.e., both are not None).

        Returns
        -------
        int
            The number of models in the best models table.

        """
        if n is not None and delta is not None:
            text = 'Cannot specify both n and delta - choose one...'
            self.logger.error(text)
            raise ValueError(text)
        elif n is None and delta is None:
            n = 10
            self.logger.info('No parameters specified - making table with '
                             '10 best models')
        if filename is None:
            path_noext, ext = os.path.splitext(self.filename)
            filename = path_noext + '_best' + ext
        else:
            filename = self.config.settings.io_settings['output_directory'] + \
                       filename
        if os.path.exists(filename):
            path_noext, ext = os.path.splitext(filename)
            backup_filename = path_noext + '_backup' + ext
            shutil.copy2(filename, backup_filename)
            self.logger.warning(f'File {filename} will be overwritten, '
                                f'backup {backup_filename} created.')
        if n is not None:
            table_best = self.get_best_n_models(which_chi2=which_chi2, n=n)
        else:
            table_best=self.get_mods_within_chi2_thresh(which_chi2=which_chi2,
                                                        delta=delta)
        table_best.write(filename, format='ascii.ecsv', overwrite=True)
        self.logger.info(f'Table of best models written to file {filename}.')
        return len(table_best)

    def remove_unused_orblibs(self):
        """
        Removes orbit libraries for 'bad' models.

        Frees disk space by deleting data, keeping only model data required
        by the ``beta_plot`` and ``mass_plot`` plotting routines.
        Keeps data of models with (kin)chi2 values less than or equal to
        ``sqrt(2 * number of kinematic observations) * min(chi2)``, but at
        least 3 models.
        Will mark a deleted orbit library with ``orblib_done=False`` and
        ``weights_done=False`` in the all_models table. If an orblib cannot
        be deleted because it is used by another model, only the nnls
        data will be deleted, which is marked by ``weights_done=False``
        in the all_models table.

        Returns
        -------
        bool
            ``True`` if data deletion has been attempted,
            ``False`` if no data to delete could be identified.

        """
        which_chi2=self.config.settings.parameter_space_settings['which_chi2']
        chi2_min = min(self.table[which_chi2])
        chi2_abs_thresh = 3 * np.sqrt(self.config.get_2n_obs())
        model_rows_keep = \
            self.get_mods_within_chi2_thresh(delta=chi2_abs_thresh)
        model_rows_del = \
            self.get_mods_within_chi2_thresh(delta=-chi2_abs_thresh)
        if len(model_rows_keep) < 3:
            self.logger.debug('Less than 3 models to keep, will keep 3 anyway.')
            model_rows_keep = self.get_best_n_models(n=3)
            model_rows_del = self.get_best_n_models(n=-3)
        self.logger.debug(f'Will remove data of {len(model_rows_del)} '
                  f'models with {which_chi2} > {chi2_min+chi2_abs_thresh}, '
                  f'keep data of {len(model_rows_keep)} models.')

        if len(model_rows_del) == 0:
            self.logger.info('Nothing to do.')
            return False

        # parameters that identify an orblib
        orblib_parameters = self.config.parspace.par_names[:]
        ml = 'ml'
        try:
            orblib_parameters.remove(ml)
        except:
            self.logger.error(f"Parameter '{ml}' not found - check "
                              "implementation")
            raise

        # now try to remove the data...
        n_removed = 0
        for model_row_del in model_rows_del:

            # get model object and row id of model whose data to delete
            parset = model_row_del[self.config.parspace.par_names]
            model = self.get_model_from_parset(parset)
            row_id = self.get_row_from_model(model)

            # remove unused orblibs
            delete_orblib = True
            for model_row_keep in model_rows_keep: # orblib used by others?
                if np.allclose(tuple(model_row_del[orblib_parameters]),
                               tuple(model_row_keep[orblib_parameters])):
                    delete_orblib = False
                    self.logger.debug("Orblib of model "
                        f"{tuple(model_row_del[orblib_parameters])} "
                        "still in use - will delete weight solving data only.")
                    break
            if delete_orblib:
                directory = model.directory_noml
                try:
                    shutil.rmtree(directory)
                    self.logger.info("Orblib of model "
                        f"{tuple(model_row_del[orblib_parameters])} "
                        f"in {directory} removed.")
                    n_removed += 1
                except:
                    self.logger.warning("Cannot remove orblib of model "
                        f"{tuple(model_row_del[orblib_parameters])} in "
                        f"{directory}, perhaps it was already removed before.")
                self.table[row_id]['orblib_done'] = False
            else:
                # orblib must be kept, but we can delete the model's nnls data
                directory = model.directory
                try:
                    shutil.rmtree(directory)
                    self.logger.info("Weight solving data of model "
                     f"{tuple(model_row_del[self.config.parspace.par_names])} "
                     f"in {directory} removed.")
                    n_removed += 1
                except:
                    self.logger.warning("Cannot remove weight solving data of "
                        f"model {tuple(model_row_del[orblib_parameters])} in "
                        f"{directory}, perhaps it was already removed before.")
            self.table[row_id]['weights_done'] = False
        self.save()
        self.logger.info(f'Removed data of {n_removed} of '
                         f'{len(model_rows_del)} identified models from disk.')
        return True


class Model(object):
    """A DYNAMITE model.

    The model can be run by running the methods (i) get_orblib, (ii) get_weights
    and (iii) (in the future) do_orbit_colouring. Running each of these methods
    will return the appropriate object, e.g. model.get_orblib() --> returns an
    OrbitLibrary object model.get_weights(...) --> returns a WeightSolver object

    Parameters
    ----------
    config : a ``dyn.config_reader.Configuration`` object
    parset : row of an Astropy Table
        contains the values of the potential parameters for this model
    directory : str
        The model directory name (without path). If None or not specified,
        the all_models_file will be searched for the directory name. If the
        all_models file does not exist, the model directory will be set to
        ``orblib_000_000/ml{ml}``.

    Returns
    -------
    Nothing returned. Attributes holding outputs are are added to the
    object when methods are run.

    """
    def __init__(self, config=None, parset=None, directory=None):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        if config is None or parset is None:
            text = f'{__class__.__name__} needs configuration object and ' \
                   'parameter set.'
            self.logger.error(text)
            raise ValueError(text)
        self.config = config
        self.check_parset(config.parspace, parset)
        self.parset = parset
        # directory of the input kinematics
        if directory is None:
           self.directory = self.get_model_directory()
        else:
           self.directory=self.config.settings.io_settings['output_directory']\
                          + 'models/' + directory
        self.logger.debug(f'Model directory string: {self.directory}')
        self.directory_noml=self.directory[:self.directory[:-1].rindex('/')+1]
        self.logger.debug('Model directory string up to ml: '
                          f'{self.directory_noml}')
        self.validate_config_file()

    def validate_config_file(self):
        """
        Validate the content of the config file against the model's config file

        Upon solving a model, DYNAMITE creates a backup of the config file in
        the model directory. Instantiating a model later (e.g., for plotting)
        using a config file that is incompatible with the one used to create
        the model can lead to problems (e.g., differently sized orbit library).
        This method validates the "global" config file against the one in the
        model directory (if existing).

        Differing config files may be ok and intended (e.g., due to expanding
        the parameter space). Therefore, the main purpose of this method is to
        add warnings to the log to alert the user.

        Raises
        ------
        FileNotFoundError
            If the "global" config file cannot be found.

        Returns
        -------
        bool
            ``False`` if a config file backup is successfully found in the
            model directory and it differs from the "global" config file.
            ``True`` otherwise (no config file backup could be identified or
            the config file backup is identical to the global config file).

        """
        if not os.path.isfile(self.config.config_file_name):
            txt = f'Unexpected: config file {self.config.config_file_name}' + \
                  f' not found (looking in {os.getcwd()}).'
            self.logger.error(txt)
            raise FileNotFoundError(txt)
        model_yaml_files = glob.glob(self.directory+'*.yaml')
        n_yaml_files = len(model_yaml_files)
        if n_yaml_files == 0:
            self.logger.debug(f'No config file backup in {self.directory} '
                              'found - probably a new model.')
            return True  # ####################
        if n_yaml_files == 1:
            f_i = 0
        else:
            try:
                f_i = model_yaml_files.index(self.directory +
                                             self.config.config_file_name)
            except ValueError:
                self.logger.warning('More than one .yaml file found in '
                                    f'{self.directory}. No file name matches '
                                    'the config file, no check possible.')
                return True  # ####################
        model_config_file_name = model_yaml_files[f_i]
        with open(self.config.config_file_name) as c_f:
            config_file = c_f.readlines()
        with open(model_config_file_name) as c_f:
            model_config_file = c_f.readlines()
        c_diff = difflib.unified_diff(config_file,
                                      model_config_file,
                                      fromfile=self.config.config_file_name,
                                      tofile=model_config_file_name,
                                      n=0)
        c_diff = list(c_diff)
        if len(c_diff) > 0:
            self.logger.warning('ACTION REQUIRED, PLEASE CHECK: '
                                'The current config file '
                                f'{self.config.config_file_name} differs from '
                                f'the config file {model_config_file_name} '
                                'backup in the model directory. Diff output:\n'
                                f'{"".join(c_diff)}.')
            return False  # ####################
        return True

    def get_model_directory(self):
        """get the name of this model's output directory

        Returns
        -------
        string

        """
        output_directory = self.config.settings.io_settings['output_directory']
        directory = output_directory + 'models/'
        models_file = output_directory \
                      + self.config.settings.io_settings['all_models_file']
        try:
            all_models = ascii.read(models_file)
            self.logger.debug(f'Setting model dir from file {models_file}...')
        except FileNotFoundError:
            sformat = self.config.system.parameters[0].sformat # ml's format
            ml_dir = f"/ml{self.parset['ml']:{sformat}}/"
            directory += f'orblib_000_000{ml_dir}'
            self.logger.info(f'The all_models file {models_file} does not '
                             f'exist - model directory set to {directory}.')
            return directory #######################################
        except:
            self.logger.error('Error reading all_models file. '
                              'Cannot set model directory.')
            raise
        for idx, parset in enumerate(all_models[self.config.parspace.par_names]):
            if np.allclose(tuple(parset),tuple(self.parset)):
                dir_string = all_models['directory'][idx]
                if dir_string == 'None': # yes, really...
                    text = f'Something went wrong reading the model # {idx} ' \
                           'directory. Model not yet computed?'
                    self.logger.error(text)
                    raise ValueError(text)
                directory += dir_string
                break
        else:
            text = f'Cannot set model directory: parset {self.parset} ' \
                   f'not found in {models_file}.'
            self.logger.error(text)
            raise ValueError(text)
        self.logger.debug(f'...model directory {directory} read from file.')
        return directory

    def create_model_directory(self, path):
        """make a directory if it does not yet exist

        Parameters
        ----------
        path : string
            directory path to make
        """
        if not os.path.exists(path):
            os.makedirs(path)
            self.logger.debug(f'Created directory {path}')
        else:
            self.logger.debug(f'Using existing directory {path}')

    def setup_directories(self):
        """setup directories
        """
        # create self.directory if it doesn't exist
        self.create_model_directory(self.directory)
        # and also the model directories
        self.create_model_directory(self.directory_noml+'infil/')
        self.create_model_directory(self.directory_noml+'datfil/')

    def get_orblib(self):
        """Make the orbit library

        Returns
        -------
        a ``dyn.orblib.OrbitLibrary`` object

        """
        orblib = dyn_orblib.LegacyOrbitLibrary(
                config=self.config,
                mod_dir=self.directory_noml,
                parset=self.parset)
        orblib.get_orblib()
        return orblib

    def get_weights(self, orblib=None):
        """Get the orbital weights

        Gets the orbital weights and chi2 values by calling the appropriate
        ``WeightSolver.solve()`` method.

        Parameters
        ----------
        orblib : a ``dyn.orblib.OrbitLibrary`` object

        Returns
        -------
        weight_solver : a ``dyn.weight_solver.WeightSolver`` object
            sets attributes:
                - ``self.weights``
                - ``self.chi2``
                - ``self.kinchi2``
                - ``self.kinmapchi2``

        """
        ws_type = self.config.settings.weight_solver_settings['type']
        if ws_type=='LegacyWeightSolver':
            weight_solver = ws.LegacyWeightSolver(
                    config=self.config,
                    directory_with_ml=self.directory)
        elif ws_type=='NNLS':
            weight_solver = ws.NNLS(
                    config=self.config,
                    directory_with_ml=self.directory)
        else:
            raise ValueError(f'Unknown WeightSolver type {ws_type}.')
        weights, chi2_tot, chi2_kin, chi2_kinmap = weight_solver.solve(orblib)
        self.chi2 = chi2_tot # instrinsic/projected mass + GH coeeficients 1-Ngh
        self.kinchi2 = chi2_kin # GH coeeficients 1-Ngh
        self.kinmapchi2 = chi2_kinmap
        self.weights = weights
        return weight_solver

    def check_parset(self, parspace, parset):
        """
        Validate a parset

        Given parameter values in parset, the validate_parspace method of
        the parameter space is executed. If a parameter exists in parspace
        but not in parset, a warning will be issued and the parameter
        will remain unchanged. If parset tries to set the value of
        a parameter not existing in parspace, an exception will be raised.
        Validating relies on exceptions raised by validate_parspace.

        Parameters
        ----------
        parspace : ``dyn.parameter_space.ParameterSpace``
            A list of parameter objects.
        parset : row of an Astropy Table
            Contains parameter values to be checked against the settings in
            parspace.

        Raises
        ------
        ValueError
            If at least one parameter in parset is unknown to parspace.

        Returns
        -------
        None.

        """
        parspace_copy = copy.deepcopy(parspace)
        parspace_par_names = [p.name for p in parspace_copy]
        unknown_pars = [p for p in parset.colnames
                          if p not in parspace_par_names]
        if len(unknown_pars) > 0:
            text = f"Parset parameters {unknown_pars} don't exist in parset."
            self.logger.error(text)
            raise ValueError(text)
        for par_idx, par_name in enumerate(parspace_par_names):
            if par_name not in parset.colnames:
                self.logger.warning(f'Parspace parameter {par_name} '
                                    'unchanged (not in parset).')
            else:
                par = parspace_copy[par_idx]
                par.par_value = parset[par_name]
        parspace_copy.validate_parspace()
        self.logger.debug('parset validated against parspace.')

# end
