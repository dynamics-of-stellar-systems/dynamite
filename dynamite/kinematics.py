import numpy as np
from scipy import special, stats
from astropy import table
import logging
import os
import h5py

from dynamite import data

class Kinematics(data.Data):
    """
    Abstract class for Kinematics

    Specific implementations (e.g. ``GaussHermite``) should be implemented as
    child classes
    """
    values = []
    def __init__(self,
                 weight=None,
                 type=None,
                 hist_width='default',
                 hist_center='default',
                 hist_bins='default',
                 **kwargs
                 ):
        super().__init__(**kwargs)
        if hasattr(self, 'data'):
            self.weight = weight
            self.type = type
            if hist_width=='default':
                self.set_default_hist_width()
            else:
                self.hist_width = hist_width
            if hist_center=='default':
                self.set_default_hist_center()
            else:
                self.hist_center = hist_center
            if hist_bins=='default':
                self.set_default_hist_bins()
            else:
                self.hist_bins = hist_bins
            self.__class__.values = list(self.__dict__.keys())
            self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
            if self.weight==None or self.type==None or self.hist_width==None or \
                    self.hist_center==None or self.hist_bins==None:
                text = 'Kinematics need (weight, type, hist_width, hist_center, '\
                   f'hist_bins), but has ({self.weight}, {self.type}, ' \
                   f'{self.hist_width}, {self.hist_center}, {self.hist_bins})'
                self.logger.error(text)
                raise ValueError(text)

    def update(self, **kwargs):
        """
        Updates one or more attributes, including consistency check

        Parameters
        ----------
        **kwargs : attribute=value pairs

        Raises
        ------
        ValueError
            If class attribute does not exist

        Returns
        -------
        None.

        """

        for k, v in kwargs.items():
            if k not in self.__class__.values:
                text = 'Invalid kinematics key ' + k + '. Allowed keys: ' + \
                    str(tuple(self.__class__.values))
                self.logger.error(text)
                raise ValueError(text)
            setattr(self, k, v)

    def validate(self): # here we can put more validation...
        if sorted(self.__class__.values) != sorted(self.__dict__.keys()):
            text = 'Kinematics attributes can only be ' + \
                str(tuple(self.__class__.values)) + ', not ' + \
                str(tuple(self.__dict__.keys()))
            self.logger.error(text)
            raise ValueError(text)

    def __repr__(self):
        return f'{self.__class__.__name__}({self.__dict__})'

    def transform_orblib_to_observables(self,
                                        losvd_histograms,
                                        weight_solver_settings):
        """trasform orbit library to observed kinematics

        This is a placeholder method. Specific implementations/subclasses of
        ``Kinematics`` should replace this with methods which transform an orbit
        library LOSVD to that particular type of kinematics.

        Parameters
        ----------
        losvd_histograms : ``dyn.kinematics.Histogram``
            the LOSVD of an orbit library
        weight_solver_settings : dict
            weight solver settings

        Returns
        -------
        object
            the orbit library LOSVD transofrmed to observed kinematics

        """
        observables = 0
        return observables

    def get_observed_values_and_uncertainties(self, weight_solver_settings):
        """Extract mean and uncertainties from the kinematic set

        This is a placeholder method. Specific implementations/subclasses of
        ``Kinematics`` should replace this with methods which extract the mean
        and uncertainties of that particular kinematic type.

        Parameters
        ----------
        weight_solver_settings : type
            Description of parameter `weight_solver_settings`.

        Returns
        -------
        tuple
            (observed_values, uncertainties)
        """
        observed_values = 0
        uncertainties = 0
        return observed_values, uncertainties

class GaussHermite(Kinematics, data.Integrated):
    """LOSVDs described by a Gauss Hermite expansion

    Using the Capellari et al 2016 convention
    """
    def __init__(self, **kwargs):
        # super goes left to right, i.e. first calls "Kinematics" __init__, then
        # calls data.Integrated's __init__
        super().__init__(**kwargs)
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        if hasattr(self, 'data'):
            self.max_gh_order = self.get_highest_order_gh_coefficient()
            self.n_apertures = len(self.data)

    def get_highest_order_gh_coefficient(self, max_gh_check=20):
        """Get max order GH coeeff from data table

        Checks the data table for columns titled ['h{i}'] and ['dh{i}'], and
        return the largest i such that both exist

        Parameters
        ----------
        max_gh_check : int, optional
            check up to this order

        Returns
        -------
        int
            the highest order GH present in the data table

        """
        colnames = self.data.colnames
        gh_order_in_table = [i for i in range(max_gh_check) if f'h{i}' in colnames and f'dh{i}' in colnames]
        try:
            max_gh_order = np.max(gh_order_in_table)
        except ValueError:
            max_gh_order = None
        return max_gh_order

    def read_file_old_format(self, filename):
        """Read the old format of GH kinematics

        Old format is that used in triaxialschwarzschild AKA schwpy codes

        Parameters
        ----------
        filename : string
            filename of old format of GH kinematics file

        Returns
        -------
        Astropy table
            the GH kinematics data

        """
        f = open(filename)
        header = f.readline()
        f.close()
        header = header.split()
        n_vbins, n_gh = header
        n_vbins, n_gh = int(n_vbins), int(n_gh)
        self.n_vbins = n_vbins
        self.n_gh = n_gh
        names = ['vbin_id', 'v', 'dv', 'sigma', 'dsigma', 'n_gh']
        for i_gh in range(3, self.n_gh+1):
            names += [f'h{i_gh}', f'dh{i_gh}']
        ncols = len(names)
        dtype = [float for i in range(ncols)]
        dtype[0] = int
        dtype[5] = int
        data = np.genfromtxt(filename,
                             skip_header=1,
                             names=names,
                             dtype=dtype)
        self.logger.debug(f'File {filename} read (old format)')
        return data

    def convert_file_from_old_format(self,
                                     filename_old_format,
                                     filename_new_format):
        """Convert old format of GH kinematics to new format

        Old format is that used in triaxialschwarzschild AKA schwpy codes

        Parameters
        ----------
        filename_old_format : string
            filename of old/shwpy format of GH kinematics file
        filename_new_format : string
            desired filename of new Astropy ECSV format

        Returns
        -------
        None
            creates a an astropy ECSV file at filename_new_format

        """
        data = self.read_file_old_format(filename_old_format)
        data = table.Table(data)
        data.remove_column('n_gh')
        data.write(filename_new_format, format='ascii.ecsv')
        self.logger.debug(f'File {filename_new_format} written (new format)')
        return

    def convert_to_old_format(self, filename_old_format):
        nbins = len(self.data)
        # extract n_gh from list of column names
        # colnames = self.data.colnames
        # idx = ['h' in x and 'd' not in x for x in colnames]
        # idx = np.where(idx)[0]
        # n_gh = [int(colnames[idx0][1]) for idx0 in idx]
        # n_gh = np.max(n_gh)
        n_gh = self.get_highest_order_gh_coefficient()
        # write comment string
        comment = '{0} {1}'.format(nbins, n_gh)
        idx = np.arange(nbins)+1
        velSym = self.data['v'].data
        dvelSym = self.data['dv'].data
        sigSym = self.data['sigma'].data
        dsigSym = self.data['dsigma'].data
        n_gh_col = np.full_like(velSym, n_gh)
        array_to_print = [idx, velSym, dvelSym, sigSym, dsigSym, n_gh_col]
        fmt = '%5i %13.13s %13.13s %13.13s %13.13s %5i '
        for i in range(3, n_gh+1):
            hi_Sym = self.data[f'h{i}'].data
            dhi_Sym = self.data[f'dh{i}'].data
            array_to_print += [hi_Sym, dhi_Sym]
            fmt += '%13.13s %13.13s '
        array_to_print = np.transpose(array_to_print)
        np.savetxt(filename_old_format,
                   array_to_print,
                   fmt = fmt,
                   header=comment,
                   comments='')
        self.logger.debug(f'File {filename_old_format} written (old format)')
        return 0

    def get_hermite_polynomial_coeffients(self, max_order=None):
        """Get Hermite poly coeffients

        Normalised as in eqn 14 of Capellari 16

        Parameters
        ----------
        max_order : int
            maximum order hermite polynomial desired e.g. max_order = 3 means
            use h0, h1, h2, h3

        Returns
        -------
        array (max_order+1, max_order+1)
            coeffients[i,j] = coef of x^j in polynomial of order i

        """
        if max_order is None:
            max_order = self.n_gh
        coeffients = []
        for i in range(0, max_order+1):
            # physicists hermite polynomials
            coef_i = special.hermite(i)
            coef_i = coef_i.coefficients
            # reverse poly1d array so that j'th entry is coeefficient of x^j
            coef_i = coef_i[::-1]
            # scale according to eqn 14 of Capellari 16
            coef_i *= (special.factorial(i) * 2**i)**-0.5
            # fill poly(i) with zeros for 0*x^j for j>i
            coef_i = np.concatenate((coef_i, np.zeros(max_order-i)))
            coeffients += [coef_i]
        coeffients = np.vstack(coeffients)
        return coeffients

    def standardise_velocities(self, v, v_mu, v_sig):
        """ Take away v_mu, divide by v_sig

        Parameters
        ----------
        v : array
            input velocity array
        v_mu : array (n_regions,)
            gauss hermite v parameters
        v_sig : array (n_regions,)
            gauss hermite sigma parameters

        Returns
        -------
        array (n_regions,) + v.shape
            velocities whitened by array v_mu, v_sigma

        """
        v = np.atleast_2d(v)
        v_mu = np.atleast_1d(v_mu)
        v_sig = np.atleast_1d(v_sig)
        assert v_mu.shape==v_mu.shape
        w = (v.T - v_mu)/v_sig
        w = w.T
        return w

    def evaluate_hermite_polynomials(self,
                                     coeffients,
                                     w,
                                     standardised=True,
                                     v_mu=None,
                                     v_sig=None):
        """Evaluate Hermite polynomials

        Parameters
        ----------
        coeffients : array (n_herm, n_herm)
            coefficients of hermite polynomials as given by method
            ``get_hermite_polynomial_coeffients``
        w : array
            if standardised==True then w is array of shape (n_regions, n_vbins)
            of standardised (AKA whitened) velocities
            else, w is array of shape (n_vbins,) of physical velocities
            and arrays v_mu and v_sig with shape (n_regions,) must be set
        standardised : Boolean
            whether or not velocities w have been standardised by aperture v/sig
        v_mu : None or array shape (n_regions,)
            aperture v_mu's
        v_sig : None or array shape (n_regions,)
            aperture v_sigma's

        Returns
        -------
        array shape (n_herm, n_regions, n_vbins)
            array[i,j,:] is the i'th Hermite polynomial evaluated at
            standardised velocities w in the j'th region

        """
        if not standardised:
            w = self.standardise_velocities(w, v_mu, v_sig)
        n_herm = coeffients.shape[0]
        w_pow_i = np.stack([w**i for i in range(n_herm)])
        # coeffients has shape (n_herm, n_herm)
        # w_pow_i has shape (n_herm, n_regions, n_vbins)
        result = np.einsum('ij,jkl->ikl', coeffients, w_pow_i, optimize=True)
        return result

    def evaluate_losvd(self, v, v_mu, v_sig, h):
        """ evaluate LOSVDs

        Evaluate the quantity

        .. math::

            \mathrm{LOSVD}(v) = \\frac{1}{v_\sigma} \mathcal{N}(w; 0, 1^2)
                \Sigma_{m=0}^{M} h_m H_m(w)

        where normalised velocity :math:`w = (v-v_\mu)/v_\sigma`

        Parameters
        ----------
        v : array
            input velocity array
        v_mu : array (n_regions,)
            gauss hermite v parameters
        v_sig : array (n_regions,)
            gauss hermite sigma parameters
        h : array (n_hists, n_regions, n_herm)
            gauss hermite expansion coefficients

        Returns
        -------
        array shape same as v
            values of gauss hermite expansion evaluated at v

        """
        w = self.standardise_velocities(v, v_mu, v_sig)
        n_herm = h.shape[2]
        max_order = n_herm - 1
        coef = self.get_hermite_polynomial_coeffients(max_order=max_order)
        nrm = stats.norm()
        hpolys = self.evaluate_hermite_polynomials(coef, w)
        losvd = np.einsum('i,ij,kil,lij->kij',
                          1./v_sig,
                          nrm.pdf(w),
                          h,
                          hpolys,
                          optimize=True)
        return losvd

    def evaluate_losvd_normalisation(self, h):
        """Evaluate LOSVD normalisation

        Evaluate the normalising integral

        .. math::

            \int_{-\infty}^{\infty} \mathrm{LOSVD}(v) dv

        which for a GH expansion is given by

        .. math::

            \Sigma_{m=0}^{M} b_m a_m

        where
            - :math:`a_m` are the coefficients of w in the polynomial :math:`\Sigma_{m=0}^{M} h_m H_m(w)`
            - :math:`b_m` =
                    - 1 if m=0
                    - 0 if m is odd
                    - (m-1)!!    if m is non-zero and even

        and !! is a 'double factorial' - which does **not** mean two factorials
        but the product integers < n with same even/odd parity as n

        Parameters
        ----------
        v_sig : array (n_regions,)
            gauss hermite sigma parameters
        h : array (n_hists, n_regions, n_herm)
            gauss hermite expansion coefficients for some number of histograms
            and regions

        Returns
        -------
        array (n_hists, n_regions)
        """
        n_herm = h.shape[2]
        max_order = n_herm - 1
        # coeffients[i,j] = coef of x^j in polynomial of order i
        coef = self.get_hermite_polynomial_coeffients(max_order=max_order)
        a_m = np.einsum('ijk,kl', h, coef)
        b_m = np.arange(0, max_order+1, 1)
        b_m = special.factorial2(b_m - 1)
        b_m[0] = 1.
        b_m[1::2] = 0.
        normalisation = np.sum(b_m * a_m, 2)
        return normalisation

    def get_gh_expansion_coefficients(self,
                                      v_mu=None,
                                      v_sig=None,
                                      vel_hist=None,
                                      max_order=4):
        """Calcuate GH expansion coeffients given an LOSVD

        Expand LOSVD around a given v_mu and v_sig using eqn 7 of
        vd Marel & Franx 93

        Parameters
        ----------
        v_mu : array (n_regions,)
            gauss hermite v parameters
        v_sig : array (n_regions,)
            gauss hermite sigma parameters
        vel_hist : Histogram object
            velocity histograms where vel_hist.y has shape
            (n_orbits, n_vbins, n_regions)
        max_order : int
            maximum order hermite polynomial desired in the expansion
            e.g. max_order = 1 --> use h0, h1
            i.e. number of hermite polys = max_order + 1

        Returns
        -------
        h : array (n_hists, n_regions, max_order+1)
            where h[i,j,k] is order k GH coeffient of histogram i in region j

        """
        assert max_order>=0
        w = self.standardise_velocities(vel_hist.x, v_mu, v_sig)
        coef = self.get_hermite_polynomial_coeffients(max_order=max_order)
        nrm = stats.norm()
        hpolys = self.evaluate_hermite_polynomials(coef, w)
        # TODO: optimize the next line for (i) vel_hist.dx is constant, (ii)
        # arrays are too large for memory e.g. using dask
        h = np.einsum('ijk,kj,lkj,j->ikl', # integral in eqn 7
                      vel_hist.y,
                      nrm.pdf(w),
                      hpolys,
                      vel_hist.dx,
                      optimize=True)
        h *= 2 * np.pi**0.5 # pre-factor in eqn 7
        return h

    def transform_orblib_to_observables(self,
                                        losvd_histograms,
                                        weight_solver_settings):
        max_gh_order = weight_solver_settings['number_GH']
        v_mu = self.data['v']
        v_sig = self.data['sigma']
        orblib_gh_coefs = self.get_gh_expansion_coefficients(
            v_mu=v_mu,
            v_sig=v_sig,
            vel_hist=losvd_histograms,
            max_order=max_gh_order)
        # in triaxnnls, GH coefficients are divided by velocity spacing of the
        # histogram. This is equivalent to normalising LOSVDs as probability
        # densities before calculating the GH coefficients. For now, do as was
        # done previously:
        dv = losvd_histograms.dx
        assert np.allclose(dv, dv[0]), 'LOSVD velocity spacing must be uniform'
        orblib_gh_coefs /= dv[0]
        # remove h0 as this is not fit
        orblib_gh_coefs = orblib_gh_coefs[:,:,1:]
        return orblib_gh_coefs

    def get_observed_values_and_uncertainties(self, weight_solver_settings):
        """Extract mean/sigma from the GH kinematics

        Parameters
        ----------
        weight_solver_settings : type
            Description of parameter `weight_solver_settings`.

        Returns
        -------
        tuple
            (observed_values, uncertainties), where:
            - observed_values is array of GH exapnsion coefficients of shape
            (n_apertures, max_gh_order)
            - uncertainties is array of uncertainties on GH exapnsion
            coefficients of shape (n_apertures, max_gh_order)

        """
        max_gh_order = weight_solver_settings['number_GH']
        # construct observed values
        observed_values = np.zeros((self.n_apertures, max_gh_order))
        # h1, h2 = 0, 0
        # h3, h4, etc... are taken from the data table
        for i in range(3, self.max_gh_order+1):
            observed_values[:,i-1] = self.data[f'h{i}']
        # construct uncertainties
        uncertainties = np.zeros_like(observed_values)
        # uncertainties on h1,h2 from vdMarel + Franx 93
        uncertainties[:,0] = self.data['dv']/np.sqrt(2)/self.data['sigma']
        uncertainties[:,1] = self.data['dsigma']/np.sqrt(2)/self.data['sigma']
        # uncertainties h3, h4, etc... are taken from data table
        for i in range(3, self.max_gh_order+1):
            uncertainties[:,i-1] = self.data[f'dh{i}']
        # add the systematic uncertainties
        systematics = weight_solver_settings['GH_sys_err']
        if type(systematics) is str:
            systematics = systematics.split(' ')
            systematics = [float(x) for x in systematics]
        systematics = np.array(systematics)
        systematics = systematics[0:max_gh_order]
        uncertainties = (uncertainties**2. + systematics**2.)**0.5
        return observed_values, uncertainties

    def set_default_hist_width(self, n_sig=3.):
        """Sets default histogram width

        Set it to

        :math:`2 * max(|v| + n_\mathrm{sig}*\sigma)`

        i.e. double the largest velcoity present in the observed LOSVD

        Parameters
        ----------
        n_sig : float
            number of sigma above mean velocity to extend the velocity histogram

        """
        v, sig = self.data['v'], self.data['sigma']
        max_abs_v_plus_3sig = np.max(np.abs(v) + n_sig*sig)
        hist_width = 2.*max_abs_v_plus_3sig
        self.hist_width = hist_width

    def set_default_hist_center(self):
        """Sets default histogram center to 0.

        """
        self.hist_center = 0.

    def set_default_hist_bins(self, f_sig=0.1):
        """Sets default nbins so they are roughly f_sig*min(sig)

        Parameters
        ----------
        f_sig : float
            fraction of minimum sigma to (approximately) set histogram bin width

        """
        sig = self.data['sigma']
        dv_approx = f_sig * np.min(sig)
        hist_bins = np.ceil(self.hist_width/dv_approx)
        hist_bins = int(hist_bins)
        # hist_bins must be odd
        if hist_bins%2==0:
            hist_bins += 1
        self.hist_bins = hist_bins

class Histogram(object):
    """LOSVD histograms

    Parameters
    ----------
    xedg : array (n_bins+1,)
        histogram bin edges
    y : (n_orbits, n_bins+1, n_apertures)
        histogram values
    normalise : bool, default=True
        whether to normalise to pdf

    Attributes
    ----------
    x : array (n_bins,)
        bin centers
    dx : array (n_bins,)
        bin widths
    normalised : bool
        whether or not has been normalised to pdf

    """
    def __init__(self, xedg=None, y=None, normalise=False):
        self.xedg = xedg
        self.x = (xedg[:-1] + xedg[1:])/2.
        self.dx = xedg[1:] - xedg[:-1]
        self.y = y
        if normalise:
            self.normalise()

    def get_normalisation(self):
        """Get the normalsition

        Calculates ``Sum_i losvd_i * dv_i``

        Returns
        -------
        float
            the normalisation

        """
        na = np.newaxis
        norm = np.sum(self.y*self.dx[na,:,na], axis=1)
        return norm

    def normalise(self):
        """normalises the LOSVDs

        Returns
        -------
        None
            resets ``self.y`` to a normalised version

        """
        norm = self.get_normalisation()
        na = np.newaxis
        tmp = self.y/norm[:,na,:]
        # where norm=0, tmp=nan. Fix this:
        idx = np.where(norm==0.)
        tmp[idx[0],:,idx[1]] = 0.
        # replace self.y with normalised y
        self.y = tmp

    def scale_x_values(self, scale_factor):
        """scale the velocity array

        scales vel array, and dv

        Parameters
        ----------
        scale_factor : float

        Returns
        -------
        None
            resets ``self.xedg``, ``self.x``, ``self.dx`` to rescaled versions

        """
        self.xedg *= scale_factor
        self.x *= scale_factor
        self.dx *= scale_factor

    def get_mean(self):
        """Get the mean velocity

        Returns
        -------
        array shape (n_orbits, n_apertures)
            mean velcoity of losvd

        """
        na = np.newaxis
        mean = np.sum(self.x[na,:,na] * self.y * self.dx[na,:,na], axis=1)
        norm = self.get_normalisation()
        mean /= norm
        return mean

    def get_sigma(self):
        """Get the velocity dispersions

        Returns
        -------
        array shape (n_orbits, n_apertures)
            velocity dispersion of losvd

        """
        na = np.newaxis
        mean = self.get_mean()
        v_minus_mu = self.x[na,:,na]-mean[:,na,:]
        var = np.sum(v_minus_mu**2. * self.y * self.dx[na,:,na],
                     axis=1)
        norm = self.get_normalisation()
        var /= norm
        sigma = var**0.5
        return sigma

class BayesLOSVD(Kinematics, data.Integrated):
    """Bayes LOSVD kinematic data

    """
    def __init__(self, **kwargs):
        # super goes left to right, i.e. first calls "Kinematics" __init__, then
        # calls data.Integrated's __init__
        super().__init__(**kwargs)
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        if hasattr(self, 'data'):
            self.convert_losvd_columns_to_one_multidimensional_column()
            self.set_mean_v_and_sig_v_per_aperture()

    def convert_losvd_columns_to_one_multidimensional_column(self):
        """Convert 1D to multi-dim columns

        ECSV files can save 1D columns, but useful to work with multi-dim
        columns for the LOSVSD. This method converts.

        Returns
        -------
        Re-sets ``self.data['losvd']`` and ``self.data['dlosvd']`` to multi-dim
        columns

        """
        nbins = self.data.meta['nbins']
        nv = self.data.meta['nvbins']
        losvd_mean = np.zeros((nbins,nv))
        losvd_sigma = np.zeros((nbins,nv))
        for j in range(nv):
            losvd_mean[:,j] = self.data[f'losvd_{j}']
            self.data.remove_column(f'losvd_{j}')
            losvd_sigma[:,j] = self.data[f'dlosvd_{j}']
            self.data.remove_column(f'dlosvd_{j}')
        self.data['losvd'] = losvd_mean
        self.data['dlosvd'] = losvd_sigma

    def convert_multidimensional_losvd_columns_to_univariate(self):
        """Convert multi-dim columns to 1D

        ECSV files can save 1D columns, but useful to work with multi-dim
        columns for the LOSVSD. This method converts.

        Returns
        -------
        Re-sets ``self.data['losvd']`` and ``self.data['dlosvd']`` to 1D columns
        called ``self.data['losvd_{i}']`` and ``self.data['dlosvd_{i}']`` for
        i = 1, ..., N_LOSVD_bins

        """
        nbins = self.data.meta['nbins']
        nv = self.data.meta['nvbins']
        losvd_mean = self.data['losvd']
        losvd_sigma = self.data['dlosvd']
        for j in range(nv):
            self.data[f'losvd_{j}'] = losvd_mean[:,j]
            self.data[f'dlosvd_{j}'] = losvd_sigma[:,j]
        self.data.remove_column('losvd')
        self.data.remove_column('dlosvd')

    def save_data_table(self, outfile=None):
        """Special save method for BayesLOSVD data.

        Handles conversion from multi-dim --> 1D columns. Should supercede the
        generic method ``dyn.data.Data.save``

        Parameters
        ----------
        outfile : string
            Name of output file

        """
        if outfile is None:
            outfile = self.datafile
            if hasattr(self, 'input_directory'):
                outfile = self.input_directory + outfile
        self.convert_multidimensional_losvd_columns_to_univariate()
        self.data.write(outfile, format='ascii.ecsv', overwrite=True)

    def load_hdf5(self, filename):
        """
        Load a hdf5 file of BAYES-LOSVD output

        Borrowed from bayes_losvd_load_hdf5.py

        Parameters
        ----------
        filename : string
            the hdf5 file of BayesLOSVD output

        Returns
        -------
        dict
            data read from BayesLOSVD output file

        """
        self.logger.info("Loading "+filename+" data")
        # Checking file exists
        if not os.path.exists(filename):
            self.logger.error("Cannot find file "+filename)
        # Open file
        self.logger.debug("# Opening file")
        f = h5py.File(filename,'r')
        # Defining output dictionary
        struct = {}
        # Filling up dictionary
        self.logger.debug("# Loading input data:")
        input_data = f['in']
        for key,values in input_data.items():
            self.logger.debug(' - '+key)
            struct[key] = np.array(values)
        if f.get("out") != None:
            self.logger.debug("# Loading Stan results:")
            output_data = f['out']
            bins_list   = list(output_data.keys())
            for idx in bins_list:
                tmp = f['out/'+idx]
                struct[int(idx)] = {}
                for key,values in tmp.items():
                    self.logger.debug(' - ['+idx+'] '+key)
                    struct[int(idx)][key] = np.array(values)
        self.logger.info("load_hdf5 is DONE")
        return struct

    def write_losvds_to_ecsv_format(self,
                                    filename=None,
                                    outfile='bayes_losvd_kins.ecsv'):
        """Convert BayesLOSVD output to ECSV

        Reads in the hdf5 BayesLOSVD output file for all spatial bins.
        Saves the median and 68% Bayesian Credible Interval of the LOSVD in each
        LOSVD bin, into an astropy ECSV file.

        Parameters
        ----------
        filename : string
            BayesLOSVD hdf5 output file for all spatial bins.
        outfile : string
            desired name of output ECSV file

        Returns
        -------
        None

        """
        result = self.load_hdf5(filename)
        dv = float(result['velscale'])
        vcent = result['xvel']
        nv = len(vcent)
        completed_bins = [i for i in result.keys() if type(i) is int]
        completed_bins = np.sort(completed_bins)
        nbins = len(completed_bins)
        losvd_mean = np.zeros((nbins,nv))
        losvd_sigma = np.zeros((nbins,nv))
        # completed_bins may have gaps i.e. some bins may not be completed, so
        # to fill arrays use a counter `i` - do not use completed_bins as index
        i = 0
        for bin in completed_bins:
            # get median LOSVD
            losvd_mean[i] = result[bin]['losvd'][2]
            # get 68% BCI
            bci_68 = result[bin]['losvd'][3] - result[bin]['losvd'][1]
            losvd_sigma[i] = 0.5*bci_68
            i += 1
        # BAYES-LOSVD returns the velocity array (and losvds) in descening order
        # let's flip them to make it easier to work with later
        vcent = vcent[::-1]
        losvd_mean = losvd_mean[:,::-1]
        losvd_sigma = losvd_sigma[:,::-1]
        # put them in a table
        data = table.Table()
        data['binID_BayesLOSVD'] = completed_bins
        # BayesLOSVD bin indexing starts at 0 and some bins may be missing, but:
        # 1) orblib_f.f90 assumes bins start at 1
        # 2) LegacyOrbitLibrary.read_orbit_base assumes that no bins are missing
        # so let's introduce a binID_dynamite which satisfies these requirements
        data['binID_dynamite'] = np.arange(nbins)+1
        for j in range(nv):
            data[f'losvd_{j}'] = losvd_mean[:,j]
            data[f'dlosvd_{j}'] = losvd_sigma[:,j]
        # add meta-data
        meta = {'dv':dv, 'vcent':list(vcent), 'nbins':nbins, 'nvbins':nv}
        data = table.Table(data, meta=meta)
        data.write(outfile, format='ascii.ecsv', overwrite=True)
        return

    def write_aperture_and_bin_files(self,
                                     filename=None,
                                     angle_deg=0.,
                                     center='max_flux',
                                     aperture_filename='aperture.dat',
                                     bin_filename='bins.dat'):
        """
        Write ``aperture.dat`` and ``bins.dat`` files

        Parameters
        ----------
        filename : string
            filename XXX_results.hdf5 of BAYES-LOSVD output (all bins combined)
        angle_deg : float
            Angle in degrees measured counter clockwise from the galaxy major
            axis to the X-axis of the input data
        center : tuple or string
            either pair of floats (x0,y0) defining center in pixel co-ordinates
            or string 'max_flux'
        aperture_filename : string
            name of aperture file
        bin_filename : string
            name of bins file

        """
        def get_pixel_info(x):
            """Get pixel information for aperture.dat given 1D image co-ords

            Given image co-ordinates x, gives (i) the likely pixel spacing dx
            assuming x are from a regular grid, and dx is the minimum spacing
            between sorted x, (ii) the minimum x in the range, (iii) the number
            of pixels, (iii) the pixel edges

            ** this will fail if there are no adjacent values in x **

            """
            ux = np.unique(x) # this is sorted
            dx = ux[1:] - ux[:-1]
            dx = np.min(dx)
            min_x, max_x = ux[0], ux[-1]
            nx = np.round((max_x - min_x)/dx) + 1
            x_cnt = min_x + np.arange(nx) * dx
            x_edg = np.concatenate((x_cnt-dx/2., [x_cnt[-1]+dx/2.]))
            x_rng = x_edg[-1] - x_edg[0]
            min_x_rng = x_edg[0]
            return min_x_rng, x_rng, x_edg, x_cnt, dx, int(nx)
        result = self.load_hdf5(filename)
        # find bins which appear in the data table
        idx = np.isin(result['binID'], self.data['binID_BayesLOSVD'])
        x = result['x'][idx]
        y = result['y'][idx]
        # center the object
        if center == 'max_flux':
            pixel_flux = result['flux'][idx]
            idx_max_flux = np.argmax(pixel_flux)
            x_center = x[idx_max_flux]
            y_center = y[idx_max_flux]
        else:
            x_center, y_center = center
        x -= x_center
        y -= y_center
        # get bin IDs
        binID_bl = result['binID'][idx]
        binID_dyn = self.map_binID_blosvd_to_binID_dynamite(binID_bl)
        # get bin centers and add to the data table
        xbin = result['xbin'] - x_center
        ybin = result['ybin'] - y_center
        self.data['xbin'] = xbin[self.data['binID_BayesLOSVD']]
        self.data['ybin'] = ybin[self.data['binID_BayesLOSVD']]
        self.save_data_table()
        # get pixel sizes
        min_x, x_rng, x_edg, x_cnt, dx, nx = get_pixel_info(x)
        min_y, y_rng, y_edg, y_cnt, dy, ny = get_pixel_info(y)
        # Create aperture.dat
        aperture_file = open(aperture_filename, 'w')
        aperture_file.write('#counter_rotation_boxed_aperturefile_version_2 \n')
        string = '\t{0:<.6f}\t{1:<.6f} \n'.format(min_x, min_y)
        aperture_file.write(string)
        string = '\t{0:<.6f}\t{1:<.6f} \n'.format(x_rng, y_rng)
        aperture_file.write(string)
        string = '\t{0:<.6f} \n'.format(angle_deg)
        aperture_file.write(string)
        string = '\t{0}\t{1} \n'.format(nx, ny)
        aperture_file.write(string)
        string = ' aperture = -(hst_pa) + 90 \n'
        aperture_file.write(string)
        aperture_file.close()
        # Write bins.dat file
        ix = np.digitize(x, x_edg)
        ix -= 1
        iy = np.digitize(y, y_edg)
        iy -= 1
        grid = np.zeros((nx, ny), dtype=int)
        grid[ix, iy] = binID_dyn
        comment_line = '#Counterrotaton_binning_version_1\n'
        grid_size = nx*ny
        first_line = '{0}\n'.format(grid_size)
        flattened_grid = grid.T.flatten()
        bins_file = open(bin_filename, 'w')
        bins_file.write(comment_line)
        bins_file.write(first_line)
        for i in range(grid_size):
            string = f'\t{flattened_grid[i]}'
            bins_file.write(string)
            # line-break every tenth line
            if i%10==9:
                bins_file.write('\n')
        bins_file.write('\n')
        bins_file.close()

    def map_binID_blosvd_to_binID_dynamite(self, binID_blosvd):
        """Map BayesLOSVD binIDs to DYNMAITE binIDs.

        Assumes that the table `self.data` has colums `binID_BayesLOSVD` and
        `binID_dynamite` which define the mapping. Any binID_blosvd with no
        corresponding binID_dynamite are given binID_dynamite=0.

        Parameters
        ----------
        binID_blosvd : array
            array of BayesLOSVD binIDs

        Returns
        -------
        type
            corresponding array of DYNMAITE binIDs

        """
        # first find binID_blosvd's which have no corresponding binID_dynamite
        idx_missing = np.isin(binID_blosvd, self.data['binID_BayesLOSVD'])
        idx_missing = np.where(idx_missing==False)
        # do the mapping - the method to do this is taken from
        # https://stackoverflow.com/q/13572448/11231128
        idx_srt_binid_blosvd = np.argsort(self.data['binID_BayesLOSVD'])
        srt_binid_blosvd = self.data['binID_BayesLOSVD'][idx_srt_binid_blosvd]
        srt_binid_dynamite = self.data['binID_dynamite'][idx_srt_binid_blosvd]
        index = np.digitize(binID_blosvd, srt_binid_blosvd, right=True)
        binID_dynamite = srt_binid_dynamite[index]
        # for missing entries, replace with 0
        binID_dynamite[idx_missing] = 0
        return binID_dynamite

    def set_default_hist_width(self, scale=2.):
        """Set orbit histogram width

        Set it to a multiple of data histogram width. Default 2 i.e. double to
        width of observed data. Sets result to attribute ``self.hist_width``

        Parameters
        ----------
        scale : float
            scale factor

        """
        vmin = self.data.meta['vcent'][0] - self.data.meta['dv']/2.
        vmax = self.data.meta['vcent'][-1] + self.data.meta['dv']/2.
        max_vabs = np.max(np.abs([vmin, vmax]))
        self.hist_width = 2. * scale * max_vabs

    def set_default_hist_center(self):
        """Sets orbit histogram center to 0
        """
        self.hist_center = 0.

    def set_default_hist_bins(self, oversampling_factor=10):
        """Set default LOSVD nbins for orblibs

        Uses velocity spacing of the data divided by oversampling_factor.
        Also forces nbins to be odd, so that central bin in 0-centered.

        Parameters
        ----------
        oversampling_factor : float
            scale factor to divide the data velocity spacing

        Returns
        ----------
        Sets the result to attribute `self.hist_bins`

        """
        data_dv = self.data.meta['dv']
        orblib_dv = self.data.meta['dv']/oversampling_factor
        orblib_nbins = self.hist_width/orblib_dv
        orblib_nbins = int(np.ceil(orblib_nbins))
        # make nbins odd so histogrammed centered on 0 (allows orbit flipping)
        if orblib_nbins % 2 == 0:
            orblib_nbins += 1
        self.hist_bins = orblib_nbins

    def set_mean_v_and_sig_v_per_aperture(self):
        """get mean and dispersion

        Returns
        -------
        creates columns ``self.data['v']`` and ``self.data['sigma']``

        """
        # the marginalised LOSVDs saved by BAYES-losvd do not sum to 1 -
        # we must account for this when calculating moments
        losvd = (self.data['losvd'].T/np.sum(self.data['losvd'], 1)).T
        v_array = self.data.meta['vcent']
        mu_v = np.sum(v_array * losvd, 1)
        var_v = np.sum((v_array - mu_v[:,np.newaxis])**2 * losvd, 1).T
        sig_v = var_v**0.5
        self.data['v'] = mu_v
        self.data['sigma'] = sig_v

    def rebin_orblib_to_observations(self, losvd_histograms):
        """Rebin orblib to velocity spacing of observations

        Creates a matrix ``f``, the fraction of the j'th orblib vbin in the i'th
        data vbin. Multiplies orblib losvd_histograms by f to rebin.

        Parameters
        ----------
        losvd_histograms : ``dyn.kinematics.Histogram``
            a histogram object of the orblib LOSVD

        Returns
        -------
        ``dyn.kinematics.Histogram``
            a orblib LOSVD re-binned to the data velocity spacing

        """
        v_cent, dv = np.array(self.data.meta['vcent']), self.data.meta['dv']
        v_edg = np.concatenate(((v_cent-dv/2.), [v_cent[-1]+dv/2.]))
        dx = losvd_histograms.dx[0]
        assert np.allclose(losvd_histograms.dx, dx), 'vbins must be uniform'
        # construct matrix to re-bin orbits to the data velocity spacing
        na = np.newaxis
        # the boudaries of the i'th data vbin:
        v_i = v_edg[:-1,na]
        v_ip1 = v_edg[1:,na]
        # the boudaries of the j'th orblib vbin:
        x_j = losvd_histograms.xedg[na,:-1]
        x_jp1 = losvd_histograms.xedg[na,1:]
        # f = the fraction of the j'th orblib vbin in the i'th data vbin:
        f1 = (x_jp1 - v_i)/dx
        f2 = (v_ip1 - x_j)/dx
        f1[f1>1] = 1
        f1[f1<0] = 0
        f2[f2>1] = 1
        f2[f2<0] = 0
        f = np.minimum(f1, f2)
        # TODO:  check if the following is faster if we use sparseness of f
        # sparse matrix multiplication won't work with einsum, but may be faster
        rebined_orbit_vel_hist = np.einsum('ijk,lj->ilk',
                                           losvd_histograms.y,
                                           f,
                                           optimize=True)
        return rebined_orbit_vel_hist

    def transform_orblib_to_observables(self,
                                        losvd_histograms,
                                        weight_solver_settings):
        losvd_histograms = self.rebin_orblib_to_observations(losvd_histograms)
        # losvd_histograms has shape (n_orbs, n_vbins, n_aperture)
        # weight solver expects (n_orbs, n_aperture, n_vbins)
        losvd_histograms = np.swapaxes(losvd_histograms, 1, 2)
        return losvd_histograms

    def get_observed_values_and_uncertainties(self, weight_solver_settings):
        """Get LOSVD mean/uncertainties

        Parameters
        ----------
        weight_solver_settings : dict

        Returns
        -------
        tuple
            (observed_values, uncertainties), where:
            - observed_values array of shape (n_aperture, n_vbins)
            - uncertainties array of shape (n_aperture, n_vbins)

        """
        # weight solver expects arrays of shape (n_aperture, n_vbins)
        observed_values = self.data['losvd'] # shape = (n_aperture, n_vbins)
        uncertainties = self.data['dlosvd'] # shape = (n_aperture, n_vbins)
        return observed_values, uncertainties

# end
