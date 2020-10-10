import data

import numpy as np
from scipy import special, stats
from astropy import table
from astropy.io import ascii

# TODO: move some of the kwargs from the init of 'Kinematics' to the init of
# higher level data classes, e.g. all Integrated objects will need aperturefile,
# binfile, maskfile, PSF.

# QUESTION: what does values list do...? Why is it outside of any method?

class Kinematics(data.Data):
    """
    Kinematics class holding attributes and methods pertaining to kinematics data
    """
    values = []
    def __init__(self,
                 weight=None,
                 type=None,
                 **kwargs
                 ):
        self.weight = weight
        self.type = type
        self.__class__.values = list(self.__dict__.keys())
        super().__init__(**kwargs)

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
                raise ValueError('Invalid kinematics key ' + k + '. Allowed keys: ' + str(tuple(self.__class__.values)))
            setattr(self, k, v)

    def validate(self): # here we can put more validation...
        if sorted(self.__class__.values) != sorted(self.__dict__.keys()):
            raise ValueError('Kinematics attributes can only be ' + str(tuple(self.__class__.values)) + ', not ' + str(tuple(self.__dict__.keys())))

    def __repr__(self):
        return f'{self.__class__.__name__}({self.__dict__})'


class GaussHermite(Kinematics, data.Integrated):

    def __init__(self, **kwargs):
        # super goes left to right, i.e. first calls "Kinematics" __init__, then
        # calls data.Integrated's __init__
        super().__init__(**kwargs)

    def read_file_old_format(self, filename):
        f = open(filename)
        header = f.readline()
        f.close()
        header = header.split()
        n_vbins, n_gh = header
        n_vbins, n_gh = int(n_vbins), int(n_gh)
        self.n_vbins = n_vbins
        self.n_gh = n_gh
        # n_gh read from the file header is currently unused and incorrect
        # TODO: fix/remove it... for now hardcode n_gh = 4
        self.n_gh = 4
        # TODO: move names/dtypes to the file header
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
        # 6th column of kin_data.dat is == (unused and incorrect?) n_gh
        # TODO: remove it
        return data

    def convert_file_from_old_format(self,
                                     filename_old_format,
                                     filename_new_format):
        data = self.read_file_old_format(filename_old_format)
        data = table.Table(data)
        data.remove_column('n_gh')
        data.write(filename_new_format, format='ascii.ecsv')
        return

    def convert_to_old_format(self, filename_old_format):
        # taken from generate_input_califa.py
        # NOTE: this only works for n_gh = 4
        # TODO: generalise for all n_gh
        nbins = len(self.data)
        n_gh = 4
        comment = '{0} {1}'.format(nbins, n_gh)
        idx = np.arange(nbins)+1
        velSym = self.data['v'].data
        dvelSym = self.data['dv'].data
        sigSym = self.data['sigma'].data
        dsigSym = self.data['dsigma'].data
        n_gh_col = np.full_like(velSym, n_gh)
        h3Sym = self.data['h3'].data
        dh3Sym = self.data['dh3'].data
        h4Sym = self.data['h4'].data
        dh4Sym = self.data['dh4'].data
        array_to_print = np.transpose([idx,
                                       velSym, dvelSym,
                                       sigSym, dsigSym,
                                       n_gh_col,
                                       h3Sym, dh3Sym,
                                       h4Sym, dh4Sym])
        fmt = '%5i %13.8s %13.8s %13.8s %13.8s %5i %13.8s %13.8s %13.8s %13.8s'
        np.savetxt(filename_old_format,
                   array_to_print,
                   fmt = fmt,
                   header=comment,
                   comments='')
        return 0

    def get_hermite_polynomial_coeffients(self, max_order=None):
        """Get coeffients for hermite polynomials normalised as in eqn 14 of
        Capellari 2016

        Parameters
        ----------
        max_order : int
            maximum order hermite polynomial desired
            e.g. max_order = 1 --> use h0, h1
            i.e. number of hermite polys = max_order + 1

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
        """

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
        """

        Parameters
        ----------
        coeffients : array (n_herm, n_herm)
            coefficients of hermite polynomials as given by method
            get_hermite_polynomial_coeffients
        w : array
            if standardised==True
                shape (n_regions, n_vbins), standardised velocities
            else
                shape (n_vbins,), physical velocities
                and arrays v_mu and v_sig with shape (n_regions,) must be set

        Returns
        -------
        array shape (n_hists, n_regions, n_vbins)
            Hermite polynomials evaluated at w in array of

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
        """ evaluate
            losvd(v) = 1/v_sig norm(w; 0, 1^2) Sum_{m=0}^{M} h_m H_m(w)
        where normalised velocity
            w = (v-v_mu)/v_sig

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
        """Evaluate the normalising integral
            int_{-inf}^{inf} losvd(v) dv
        which is given by
            Sum_{m=0}^{M} b_m a_m
        where:
        - a_m are the coefficients of w in the polynomial
                Sum_{m=0}^{M} h_m H_m(w)
        -       { 1          if m=0
          b_m = { 0          if m is odd
                { (m-1)!!    if m is non-zero and even
        and !! is a 'double factorial' - which does *NOT* mean two factorials
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
        """Calcuate coeffients of gauss hermite expansion of histogrammed LOSVD
        around a given v_mu and v_sig i.e. evaluate qn 7 of vd Marel & Franx 93

        Parameters
        ----------
        v_mu : array (n_regions,)
            gauss hermite v parameters
        v_sig : array (n_regions,)
            gauss hermite sigma parameters
        vel_hist : Histogram object
            velocity histograms where vel_hist.y has shape (n_hists, n_vbins)
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
        h = np.einsum('ij,kj,lkj,j->ikl',           # integral in eqn 7
                      np.atleast_2d(vel_hist.y),
                      nrm.pdf(w),
                      hpolys,
                      vel_hist.dx,
                      optimize=True)
        h *= 2 * np.pi**0.5                         # pre-factor in eqn 7
        return h

    def transform_orbits_to_observables(self, orb_lib):
        # actual code to transform orbits to GH coefficients
        return gauss_hermite_coefficients


class Histogram(object):
    """Class to hold histograms

    Parameters
    ----------
    xedg : array (n_bins+1,)
        histogram bin edges
    y : (n_histograms, n_bins+1,)
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
    def __init__(self, xedg=None, y=None, normalise=True):
        self.xedg = xedg
        self.x = (xedg[:-1] + xedg[1:])/2.
        self.dx = xedg[1:] - xedg[:-1]
        self.y = y
        if normalise:
            self.normalise()
        else:
            self.normalised = False

    def normalise(self):
        norm = np.sum(self.y*self.dx, axis=-1)
        self.y = (self.y.T/norm).T
        self.normalised = True



# end
