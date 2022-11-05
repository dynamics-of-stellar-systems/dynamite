#import logging
import numpy as np
import pyfort_GaussHerm
from astropy.io import ascii
import subprocess
import os
import matplotlib.pyplot as plt
#from astropy.io import fits
#import glob
from lmfit import minimize, Parameters
from scipy.io import FortranFile
import dynamite as dyn

from plotbin.display_pixels import display_pixels
#from matplotlib.patches import Ellipse
#from scipy.optimize import minimize_scalar
#import plotbin.display_pixels as dp
#import plotbin.display_bins_generators as dbg
#from scipy import special, stats


def scale_factor_find(ml_dir):
        """Get ``ml`` of original orblib with shared parameters

        The original ``ml`` is required to rescale orbit libraries for rescaled
        potentials. This method reads it from the first entry of 9th from bottom
        line of the file ``infil/parameters_pot.in``

        Returns
        -------
        float
            the original ``ml``

        """
        infile = ml_dir + '/nn.in'
        lines = [line.rstrip('\n').split() for line in open(infile)]
        scale_factor = float((lines[-4])[0])
        return scale_factor




####################################################
def gaussfunc_gh(paramsin,x):
    amp=paramsin['amp'].value
    center=paramsin['center'].value
    sig=paramsin['sig'].value
    c1=-np.sqrt(3);
    c2=-np.sqrt(6)
    c3=2/np.sqrt(3);
    c4=np.sqrt(6)/3;
    c5=np.sqrt(6)/4
    skew=paramsin['skew'].value
    kurt=paramsin['kurt'].value
    g=(x-center)/sig
    gaustot_gh=amp*np.exp(-.5*g**2)*(1+skew*(c1*g+c3*g**3)+ kurt*(c5+c2*g**2+c4*(g**4)))

    return gaustot_gh




def read_orbit_base(wdir, fileroot, n_apertures):
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
        os.chdir(wdir)
        # unzip orblib to a temproary file with ml value attached
        # ml value is needed to prevent different processes clashing

        aaa = subprocess.run(['bunzip2', '-c', f'datfil/{fileroot}.dat.bz2'],
                             stdout=subprocess.PIPE)
        tmpfname = f'datfil/{fileroot}.dat'
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
        nconstr = tmp[0][0]
        nvhist = tmp[1][0]
        dvhist = tmp[2][0]
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
            for k in range(nconstr):
                ivmin, ivmax = orblibf.read_ints(np.int32)
                if ivmin <= ivmax:
                    tmp = orblibf.read_reals(float)
                    velhist[j, ivmin+nvhist:ivmax+nvhist+1, k] = tmp
        orblibf.close()
        os.chdir(cur_dir)
        vedg_pos = np.arange(1, nbins_vhist+1, 2) * dvhist/2.
        vedg_neg = -vedg_pos[::-1]
        vedg = np.concatenate((vedg_neg, vedg_pos))
        velhist = dyn.kinematics.Histogram(xedg=vedg,
                                           y=velhist,
                                           normalise=False)
        return velhist, density_3D




def duplicate_flip_and_interlace_orblib(orblib):
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

        error_msg = 'velocity array must be symmetric'
        assert np.allclose(orblib.xedg, -orblib.xedg[::-1]), error_msg

        losvd = orblib.y
        n_orbs, n_vel_bins, n_spatial_bins = losvd.shape
        reveresed_losvd = losvd[:, ::-1, :]
        new_losvd = np.zeros((2*n_orbs, n_vel_bins, n_spatial_bins))
        new_losvd[0::2] = losvd
        new_losvd[1::2, :] = reveresed_losvd
        new_orblib = dyn.kinematics.Histogram(xedg=orblib.xedg,
                                              y=new_losvd,
                                              normalise=False)
        return new_orblib




def combine_orblibs( orblib1, orblib2):
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

        error_msg = 'orblibs have different number of velocity bins'
        assert n_vel_bins1==n_vel_bins2, error_msg

        error_msg = 'orblibs have different velocity arrays'
        assert np.array_equal(orblib1.x, orblib2.x), error_msg
        error_msg = 'orblibs have different number of spatial bins'
        assert n_spatial_bins1==n_spatial_bins2, error_msg

        new_losvd = np.zeros((n_orbs1 + n_orbs2,
                              n_vel_bins1,
                              n_spatial_bins1))
        new_losvd[:n_orbs1] = orblib1.y
        new_losvd[n_orbs1:] = orblib2.y
        new_orblib = dyn.kinematics.Histogram(xedg=orblib1.xedg,
                                              y=new_losvd,
                                              normalise=False)
        return new_orblib




def read_losvd_histograms(wdir, n_apertures,scale_factor):
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
        print('Reading tube orbits')
        tube_orblib, tube_density_3D = read_orbit_base(wdir, 'orblib', n_apertures)
        # tube orbits are mirrored/flipped and used twice
        print('Duplicating tube orbits')
        tube_orblib = duplicate_flip_and_interlace_orblib(tube_orblib)
        # duplicate and interlace tube_density_3D array
        tube_density_3D = np.repeat(tube_density_3D, 2, axis=0)
        # read box orbits
        print('Reading box orbits')
        box_orblib, box_density_3D = read_orbit_base(wdir, 'orblibbox', n_apertures)
        # combine orblibs
        print('Combining orbits')
        orblib = combine_orblibs(tube_orblib, box_orblib)
        # combine density_3D arrays
        density_3D = np.vstack((tube_density_3D, box_density_3D))



        scale_factor = scale_factor
        print('Scaling hists')
        orblib.scale_x_values(scale_factor)
        losvd_histograms = orblib
        intrinsic_masses = density_3D
        n_orbs = losvd_histograms.y.shape[0]
        projected_masses = np.sum(losvd_histograms.y, 1)
        return losvd_histograms, intrinsic_masses, n_orbs, projected_masses


def read_free_file(input_file, n_head_element, n_row, n_column):
    arr_out=[]
    with open(input_file) as f:
        for line in f:
            for x in line.split():
                arr_out.append(float(x))
    del arr_out[0:n_head_element]
    #print(len(arr_out))
    arr_out=np.reshape(arr_out, (n_row, n_column))
    #print(arr_out)
    #print(arr_out.shape)
    return arr_out



def conv_gh_fit(vhist, b):
    p_gh=Parameters()
    p_gh.add('amp',value=np.max(vhist),vary=True);
    p_gh.add('center',value=np.percentile(b, 50),min=np.min(b),max=np.max(b));
    #p_gh.add('sig',value=(b[w]-np.min(b))/2,min=b[w]-np.min(b),max=np.max(b)-b[w]);
    p_gh.add('sig',value=np.abs(np.percentile(b,68)-np.percentile(b, 50)),min=0,max=700);
    p_gh.add('skew',value=0,vary=True,min=None,max=None);
    p_gh.add('kurt',value=0,vary=True,min=None,max=None);
    gausserr_gh = lambda p,x,y: gaussfunc_gh(p,x)-y
    fitout_gh=minimize(gausserr_gh,p_gh,args=(b,vhist))

    fitted_p_gh = fitout_gh.params
    pars_gh=[fitout_gh.params['amp'].value,
            fitout_gh.params['center'].value,
            fitout_gh.params['sig'].value,
            fitout_gh.params['skew'].value,
            fitout_gh.params['kurt'].value]
    fit_gh=gaussfunc_gh(fitted_p_gh,b)
    pars_gh=pars_gh
    resid_gh=fit_gh-vhist
    return pars_gh[1], pars_gh[2]



def conv_moments(vhis, b, dv_obs):
    L_obs_nn = a=np.where(vhis > 0, vhis, 0.0)

    mu0=np.sum(L_obs_nn, dtype=np.double)*dv_obs
    mu1=np.sum(b*L_obs_nn, dtype=np.double)*dv_obs
    mu2=np.sum((b**2)*L_obs_nn, dtype=np.double)*dv_obs
    mu3=np.sum((b**3)*L_obs_nn, dtype=np.double)*dv_obs
    mu4=np.sum((b**4)*L_obs_nn, dtype=np.double)*dv_obs

    v = mu1 / mu0
    s = np.sqrt(mu0 * mu2 - mu1 ** 2) / mu0
    return v, s



def conv_fortran(vhis,histbinsize,wbin):
    gamm, vmi, sgi, chi2 = pyfort_GaussHerm.gaussfit(veltemp=vhis,
            nvhist=wbin, nvmax=wbin, dvhist=histbinsize, nmc=10)
    HH = pyfort_GaussHerm.getgauher(vm=vmi, sg=sgi, veltemp=vhis,
            nvhist=wbin, nvmax=wbin, dvhist=histbinsize, nhermmax=4,
            gam=1.0, ingam=1)

    return vmi, sgi



def conv_losvd(orblib, nkin,xedg):

    mod_losvd = dyn.kinematics.Histogram(y=orblib[np.newaxis,:,:], xedg=xedg)
    mod_v = mod_losvd.get_mean()
    mod_sig = mod_losvd.get_sigma()

    return mod_v[0, nkin], mod_sig[0, nkin]



def comps_aphist(w_dir, losvd_hist, comps, gal, conversion):
    wdir=w_dir+str(gal)+'/figure_nn/'
    n_orbs, bins, k_bins=losvd_hist.y.shape
    #w = int(head1[0]) #bins
    #n = int(head1[1]) #apertures
    histbinsize =losvd_hist.dx
    print('shape histbin',histbinsize.shape)

    b = (np.arange(bins , dtype=float) - ((bins+1)/2)) * histbinsize

    comps=['disk', 'thin_d', 'warm_d', 'bulge', 'all']
    xedg=losvd_hist.xedg
    losvd_orbs=losvd_hist.y
    print(losvd_orbs.shape)
    for k in range(len(comps)):
        file=wdir+comps[k]+'_orb_s22.out'
        norbout, ener, i2, i3, regul, orbtype, orbw, lcut, ntot \
            = dyn.plotter.Plotter().readorbout(file)
        norm_w=np.sum(orbw)
        #print('tot orb_weight', norm_w)
        orbw=orbw
        orb_sel=norbout
        #print('selected orbs',len(orb_sel))
        losvd=np.copy(losvd_orbs)
        #print('before zeroing some orbits', np.sum(np.sum(losvd,0)))
        check=0
        #if k==3: print('orbs', orb_sel)
        for orb in range(n_orbs):

            if orb in orb_sel:
                #print('orb #', orb)
                lk=np.where(norbout==orb)[0][0]
                for i in range(losvd.shape[1]):
                    losvd[orb,i,:]=losvd[orb,i,:]*orbw[lk]
                #if check==1:print('check orbit taken',losvd[orb,:,0],losvd[orb,:,1])
            else:
                for i in range(losvd.shape[1]):

                    losvd[orb,i,:]=0.0
            check+=1

        losvd=np.sum(losvd, 0) #sum the losvd in all the orbits
        v_calculate = np.zeros(k_bins, dtype=float)
        s_calculate = np.zeros(k_bins, dtype=float)
        lsb = np.zeros(k_bins, dtype=float)

        #not calculating h3 and h4 at the moment



        for i in range(0, k_bins):
            vhis = losvd[:, i] #losvd for each kinematic bin
            dv_obs= (np.max(b)-np.min(b))/np.double(len(b))

            if conversion=='gh_fit':
                #other options are 'losvd_vsig', 'fortran', 'moments'
                v_i, sigma_i= conv_gh_fit(vhis, b)

            if conversion=='losvd_vsig':

                v_i, sigma_i= conv_losvd(losvd, i,xedg)

            if conversion=='fortran':
                v_i, sigma_i=conv_fortran(vhis, histbinsize, (bins-1)/2)

            if conversion=='moments':
                v_i, sigma_i=conv_moments(vhis, b, dv_obs)

            L_obs_nn = a=np.where(vhis > 0, vhis, 0.0)
            mu0=np.sum(L_obs_nn, dtype=np.double)*dv_obs

            lsb[i] = mu0
            v_calculate[i] = v_i
            s_calculate[i] = sigma_i

        file_o=wdir+comps[k]+'_vs_kinem_'+conversion+'.out'

        with open(file_o, 'w') as outfile:
            outfile.write(("%13.3f" % k_bins) + ("%8i" % 2 ) + '\n')
            for i in range(0, k_bins):
                if i==0:
                    print(comps[k],lsb[i],v_calculate[i],s_calculate[i])
                outfile.write(("%12i" % (i+1))  + ("%15.7f" % lsb[i]) +  ("%15.7f" % lsb[i]) +
                              ("%15.7f" % v_calculate[i]) + ("%15.7f" % v_calculate[i])  +
                              ("%15.7f" % 5.0)  + ("%15.7f" % s_calculate[i]) + ("%15.7f" % s_calculate[i]) +
                              ("%15.7f" % 5.0)  + ("%15.7f" % 0.0) + ("%15.7f" % 0.0)+ ("%15.7f" % 0.1) +
                              ("%15.7f" % 0.0)  + ("%15.7f" % 0.0) + ("%15.7f" % 0.1) + '\n')





def create_orbital_component_files_giu(w_dir=None,object=None, rootname=None, diri_o=None, ml_str = None,
                   ocut=None, Rmax_arcs=None, xrange=None):

    plotter = dyn.plotter.Plotter()
    if not ocut:
        ocut = [0.8, 0.25, -0.25]  #selection in lambda_z following Santucci+22
    print("cut lines are:", ocut)

    nn = ml_str + '/nn' + rootname
    wdir = w_dir + object + '/' + diri_o + '/' #bestfit folder
    savedir=w_dir+object+'/figure_nn/' #where to save the diles

    if not xrange:
        xrange = [0.0, Rmax_arcs]

    dir0  = wdir+ nn
    file4 = wdir+ nn + '_orb.out'
    file2 = wdir+ 'datfil/orblib.dat_orbclass.out'  #orbitlibraries
    file3 = wdir+ 'datfil/orblibbox.dat_orbclass.out'
    file3_test = os.path.isfile(file3)
    if not file3_test: file3= '%s' % file2

    mgepar, distance, th_view, ph_view, psi_view, ml, bhmass,softlen, \
        nre, lrmin, lrmax, nrth, nrrad, ndither, vv1_1, vv1_2,dm1,dm2, \
        conversion_factor,grav_const_km,parsec_km, rho_crit \
        = plotter.triaxreadparameters(w_dir=wdir)

    norb = int(nre * nrth * nrrad)

    nrow = norb
    ncol = ndither ** 3
    #print('norb', norb)
    orbclass1 = plotter.readorbclass(file=file2, nrow=norb, ncol=ndither ** 3)
    orbclass2 = plotter.readorbclass(file=file3, nrow=norb, ncol=ndither ** 3)

    #print('norb, ndither', np.max(norb), ndither)
    norbout, ener, i2, i3, regul, orbtype, orbw, lcut, ntot \
        = plotter.readorbout(filename=file4)

    #print('ener, i2, i3', np.max(ener), np.max(i2), np.max(i3))
    #print('Maxmin and Minimum(ener)', np.max(ener), np.min(ener))

    orbclass = np.dstack((orbclass1, orbclass1, orbclass2))
    print(len(orbclass))
    orbclass1a = np.copy(orbclass1)
    orbclass1a[0:3, :, :] *= -1  # the reverse rotating orbits of orbclass

    for i in range(int(0), norb):
        orbclass[:, :, i * 2] = orbclass1[:, :, i]
        orbclass[:, :, i * 2 + 1] = orbclass1a[:, :, i]

    ## define circularity of each orbit [nditcher^3, norb]
    lz = (orbclass[2, :, :] / orbclass[3, :, :] / np.sqrt(orbclass[4, :, :]))  # lambda_z = lz/(r * Vrms)
    lx = (orbclass[0, :, :] / orbclass[3, :, :] / np.sqrt(orbclass[4, :, :]))  # lambda_x = lx/(r * Vrms)
    l = (np.sqrt(np.sum(orbclass[0:3, :, :] ** 2, axis=0)) / orbclass[3, :, :] / np.sqrt(orbclass[4, :, :]))
    r = (orbclass[3, :, :] / conversion_factor)  # from km to kpc

    # average values for the orbits in the same bundle (ndither^3).
    # Only include the orbits within Rmax_arcs

    rm = np.sum(orbclass[3, :, :] / conversion_factor, axis=0) / ndither ** 3
    lzm = np.sum(np.abs(lz), axis=0) / ndither ** 3
    lxm=np.sum(lx, axis=0) / ndither ** 3
    #print("check 1", lzm, lxm)
    #s = np.ravel(np.where((rm > xrange[0]) & (rm < xrange[1])))

    # flip the sign of lz to confirm total(lz) > 0
    t = np.ravel(np.argsort(rm))

    yy = np.max(np.ravel(np.where(np.cumsum(orbw[t]) <= 0.5)))
    k = t[0:yy]
    if np.sum(np.sum(lz[:, k], axis=0) / (ndither ** 3) * orbw[k]) < 0: lz *= -1.0


    lzm_sign= np.sum(lz, axis=0) / ndither ** 3
    lxm_sign= np.sum(lx, axis=0) / ndither ** 3
    #print("check 2 - sign", lzm_sign, lxm_sign)

    print('*******Creat nn_orb.out for different bulk of orbits***********')

    norbout, ener, i2, i3, regul, orbtype, orbw, lcut, ntot \
        = plotter.readorbout(filename=file4)

    print('#orbs', len(norbout))
    ### COLD COMPONENT
    hlz=np.ravel(np.where((lzm_sign) >= ocut[0]))
    print("check 3 - cold comp", len( hlz) )
    #print("and len(ocut)", len(ocut))

    if len(ocut) == 1:
        hlz=np.ravel(np.where((lzm_sign) >= ocut[0]))
    orbw_hlz = np.copy(orbw)

    nohlz= np.ravel(np.where((lzm_sign) < ocut[0]))
    if len(ocut) == 1:
        nohlz=np.ravel(np.where((lzm_sign) < ocut[0]))
    orbw_hlz[nohlz] = -999

    with open(savedir + 'thin_d' + '_orb_s22.out', 'w') as outfile:
        norbw = len(orbw)
        outfile.write(("%12i" % len(hlz)) + '\n')
        for i in range(0, norbw):
            if orbw_hlz[i]!=-999:
                outfile.write(("%8i" % i)  + ("%8i" % ener[i]) +  ("%8i" % i2[i]) +
                          ("%8i" % i3[i]) + ("%8i" % regul[i])  +
                          ("%8i" % orbtype[i])  + ("%13.5e" % orbw_hlz[i]) +
                              ("%13.5e" % lcut[i]) +  '\n')

    ### WARM COMPONENT
    wlz= np.ravel(np.where((((lzm_sign) > ocut[1])) & ((lzm_sign)< ocut[0])))
    orbw_wlz = np.copy(orbw)
    nowlz= np.ravel(np.where(((lzm_sign) > ocut[0]) | ((lzm_sign) <= ocut[1])))
    orbw_wlz[nowlz] = -999
    print("check 3 - warm comp", len( wlz) )
    with open(savedir + 'warm_d' + '_orb_s22.out', 'w') as outfile:
        norbw = len(orbw)
        outfile.write(("%12i" % len(wlz)) + '\n')
        for i in range(0, norbw):
            if orbw_wlz[i]!=-999:
                outfile.write(("%8i" % i)  + ("%8i" % ener[i]) +  ("%8i" % i2[i]) +
                          ("%8i" % i3[i]) + ("%8i" % regul[i])  +
                          ("%8i" % orbtype[i])  + ("%13.5e" % orbw_wlz[i]) +
                              ("%13.5e" % lcut[i]) +'\n')


    ### HOT_noCR COMPONENT
    zlz= np.ravel(np.where((((lzm_sign) > ocut[2])) & ((lzm_sign)< ocut[1])))

    orbw_zlz = np.copy(orbw)

    nozlz= np.ravel(np.where(((lzm_sign) > ocut[1]) | ((lzm_sign) <= ocut[2])))
    orbw_zlz[nozlz] = -999
    print("check 3 - hot comp", len( zlz) )
    with open(savedir + 'bulge' + '_orb_s22.out', 'w') as outfile:
        norbw = len(orbw)
        outfile.write(("%12i" % len(zlz)) + '\n')
        for i in range(0, norbw):
            if orbw_zlz[i]!=-999:
                outfile.write(("%8i" % i)  + ("%8i" % ener[i]) +  ("%8i" % i2[i]) +
                          ("%8i" % i3[i]) + ("%8i" % regul[i])  +
                          ("%8i" % orbtype[i])  + ("%13.5e" % orbw_zlz[i]) +
                              ("%13.5e" % lcut[i]) +'\n')


 ### DISK COMPONENT
    dlz= np.ravel(np.where(((lzm_sign)> ocut[1])))
    orbw_dlz = np.copy(orbw)

    nodlz= np.ravel(np.where(((lzm_sign) <= ocut[1]) ))
    orbw_dlz[nodlz] = -999
    print("check 3 - disk comp", len( dlz) )

    with open(savedir + 'disk' + '_orb_s22.out', 'w') as outfile:
        norbw = len(orbw)
        outfile.write(("%12i" % len(dlz)) + '\n')
        for i in range(0, norbw):
            if orbw_dlz[i]!=-999:
                outfile.write(("%8i" % i)  + ("%8i" % ener[i]) +  ("%8i" % i2[i]) +
                          ("%8i" % i3[i]) + ("%8i" % regul[i])  +
                          ("%8i" % orbtype[i])  + ("%13.5e" % orbw_dlz[i]) +
                              ("%13.5e" % lcut[i]) +'\n')

  ### Whole COMPONENT

    orbw_allz = np.copy(orbw)

    print("check 3 - whole comp", len(orbw_allz) )
    with open(savedir + 'all' + '_orb_s22.out', 'w') as outfile:
        norbw = len(orbw)
        outfile.write(("%12i" % len(orbw)) + '\n')
        for i in range(0, norbw):
            if orbw_allz[i]!=-999:
                outfile.write(("%8i" % i)  + ("%8i" % ener[i]) +  ("%8i" % i2[i]) +
                          ("%8i" % i3[i]) + ("%8i" % regul[i])  +
                          ("%8i" % orbtype[i])  + ("%13.5e" % orbw_allz[i]) +
                              ("%13.5e" % lcut[i]) +'\n')

    tot_orb=np.sum(orbw_hlz[orbw_hlz!=-999])+ np.sum(orbw_wlz[orbw_wlz!=-999])+np.sum(orbw_zlz[orbw_zlz!=-999])
    tot_orb2=np.sum(orbw_zlz[orbw_zlz!=-999])+ np.sum(orbw_dlz[orbw_dlz!=-999])
    tot_orb_tot=np.sum(orbw_allz)
    #print('*************************')
    #print('total flux 1 and 2 (disk+bulge), and total', tot_orb, tot_orb2, tot_orb_tot)
    #print('*************************')




    #print('*************************')
    #print('total weights of hlz (thin), wlz (thick), zlz (bulge), dlz (disk), whole components are:',
    #      np.sum(orbw_hlz[orbw_hlz!=-999])/tot_orb, np.sum(orbw_wlz[orbw_wlz!=-999])/tot_orb,np.sum(orbw_zlz[orbw_zlz!=-999])/tot_orb, np.sum(orbw_dlz[orbw_dlz!=-999])/tot_orb, np.sum(orbw_allz)/tot_orb)
    #print('************************')
    #print('*************************')
    #print('total weights of hlz (cold), wlz (warm), zlz (hot), clz (CC) components are:',
    #     np.sum(orbw_hlz), np.sum(orbw_wlz),np.sum(orbw_zlz), np.sum(orbw_clz))
    #print('************************')

    return 1




#def plot_comps_giu(w_dir=None, gal=None, rootname=None
#                        ,savedata=True, xlim=None , ylim=None, Re=None, conversion=conversion,ext=""):
def plot_comps_giu(conversion, w_dir=None, gal=None, rootname=None
                        ,savedata=True, xlim=None , ylim=None, Re=None, ext=""):


    ap_file=gal+'_aperture.dat'
    bin_file=gal+'_bins.dat'

    wdir = w_dir + gal + '/figure_nn/'
    wwdir=w_dir +  gal+'/'




    #READING KINEMATIC DATA GENERATED BY triax_makevs_from_aphist
    ##############################################################
    ## COLD COMPONENT
    kinem_file = wdir + 'thin_d_vs_kinem_'+conversion+'.out'
    weight_file = wdir + 'thin_d_orb_s22.out'
    kinem_matrix = np.genfromtxt(kinem_file, skip_header=1)
    weight_matrix = np.genfromtxt(weight_file, skip_header=1)
    id, flux, flux_thin, vel_thin, veld, dvel, sig_thin, sigd, dsig, h3_thin, h3d, dh3, h4_thin, h4d, dh4 = kinem_matrix.T

    n,n1,n2,n3,n0,nc,wthin,lcutw = weight_matrix.T


    ## WARM COMPONENT
    kinem_file = wdir + 'warm_d_vs_kinem_'+conversion+'.out'
    weight_file = wdir + 'warm_d_orb_s22.out'
    kinem_matrix = np.genfromtxt(kinem_file, skip_header=1)
    weight_matrix = np.genfromtxt(weight_file, skip_header=1)
    id, flux, flux_thick, vel_thick, veld, dvel, sig_thick, sigd, dsig, h3_thick, h3d, dh3, h4_thick, h4d, dh4 = kinem_matrix.T
    n,n1,n2,n3,n0,nc,wthick,lcutthi = weight_matrix.T


   ### CC COMPONENT
    kinem_file = wdir + 'disk_vs_kinem_'+conversion+'.out'
    weight_file = wdir + 'disk_orb_s22.out'

#
    kinem_matrix = np.genfromtxt(kinem_file, skip_header=1)
    weight_matrix = np.genfromtxt(weight_file, skip_header=1)
    id, flux, flux_disk, vel_disk, veld, dvel, sig_disk, sigd, dsig, h3_disk, h3d, dh3, h4_disk, h4d, dh4 = kinem_matrix.T
    n,n1,n2,n3,n0,nc,wdisk,lcutdisk = weight_matrix.T
    #print('disk flux', np.sum(flux_disk))

    ## HOT_cr COMPONENT
    kinem_file = wdir + 'bulge_vs_kinem_'+conversion+'.out'
    weight_file = wdir + 'bulge_orb_s22.out'
    kinem_matrix = np.genfromtxt(kinem_file, skip_header=1)
    weight_matrix = np.genfromtxt(weight_file, skip_header=1)
    id, flux, flux_bulge, vel_bulge, veld, dvel, sig_bulge, sigd, dsig, h3_bulge, h3d, dh3, h4_bulge, h4d, dh4 = kinem_matrix.T
    n,n1,n2,n3,n0,nc,wbulge,lcutbul = weight_matrix.T
    #print('bulge flux', np.sum(flux_bulge))

    ###WHOLE component
    kinem_file = wdir + 'all_vs_kinem_'+conversion+'.out'
    weight_file = wdir + 'all_orb_s22.out'
    kinem_matrix = np.genfromtxt(kinem_file, skip_header=1)
    weight_matrix = np.genfromtxt(weight_file, skip_header=1)
    id, flux, flux_all, vel_all, veld, dvel, sig_all, sigd, dsig, h3_all, h3d, dh3, h4_all, h4d, dh4 = kinem_matrix.T
    n,n1,n2,n3,n0,nc,wall,lcutal = weight_matrix.T

    ##############################################################
    # Read aperture1.dat
    # The angle that is saved in this file is measured counter clock-wise
    # from the galaxy major axis to the X-axis of the input data.
    lines = [line.rstrip('\n').split() for line in open(wwdir + ap_file)]
    strhead = lines[0]
    minx =float(lines[1][0]); miny = float(lines[1][1])
    sx = float(lines[2][0]); sy = float(lines[2][1])
    #minx = -24.5; miny= -24.5#np.float(lines[1][0]); miny = np.float(lines[1][1])
    #sx=50;sy=50#sx = np.float(lines[2][0]); sy = np.float(lines[2][1])
    maxx = sx + minx; sy = sy + miny
    angle_deg = float(lines[3][0])  # - 90.0 + 180
    nx = np.int(lines[4][0]); ny = np.int(lines[4][1])

    dx = sx / nx

    print("Pixel grid dimension is dx,nx,ny,angle , ", dx, nx, ny, angle_deg)
    grid = np.zeros((nx, ny), dtype=int)

    xr = (np.arange(nx, dtype=float) * dx + minx + 0.5 * dx)
    yc = (np.arange(ny, dtype=float) * dx + miny + 0.5 * dx)

    xi = np.outer(xr, (yc * 0 + 1)); xt = xi.T.flatten()
    yi = np.outer((xr * 0 + 1), yc); yt = yi.T.flatten()

    radeg = 57.2958
    #print('PA: ', angle_deg)
    xi=xt
    yi=yt
    ##############################################################
    # Read bins1.dat
    lines_bins = [line.rstrip('\n').split() for line in open(wwdir + bin_file)]
    i = 0 ; str_head = []
    i_var = [] ; grid = []
    while i < len(lines_bins):
        for x in lines_bins[i]:
            if i == 0:
                str_head.append(str(x))
            if i == 1:
                i_var.append(np.int(x))
            if i > 1:
                grid.append(np.int(x))
        i += 1
    str_head = str(str_head[0])
    i_var = int(i_var[0])
    grid = np.ravel(np.array(grid))
    # bins start counting at 1 in fortran and at 0 in idl:
    grid = grid - 1
    ##############################################################

    s = np.ravel(np.where((grid >= 0) & (np.abs(xi) <= xlim) & (np.abs(yi) <= ylim)))
    s_wide = np.ravel(np.where(grid >= 0))
    #normalise flux
    fhist, fbinedge = np.histogram(grid[s_wide], bins=len(flux_thin))
    flux_thin = flux_thin / fhist
    fhist, fbinedge = np.histogram(grid[s_wide], bins=len(flux_thick))
    flux_thick = flux_thick / fhist

    fhist, fbinedge = np.histogram(grid[s_wide], bins=len(flux_disk))
    flux_disk = flux_disk / fhist


    fhist, fbinedge = np.histogram(grid[s_wide], bins=len(flux_bulge))
    flux_bulge = flux_bulge / fhist

    fhist, fbinedge = np.histogram(grid[s_wide], bins=len(flux_all))
    flux_all= flux_all / fhist



    tx = np.zeros((nx, ny), dtype=float)
    ty = np.zeros((nx, ny), dtype=float)


    tthin  = flux_thin[grid]
    tthick = flux_thick[grid]

    tdisk  = flux_disk[grid]
    tbulge =flux_bulge[grid]
    tall =flux_all[grid]

    #print('before normalization th  tw  tz tc thc tcw are:', np.sum(th),  np.sum(tw), np.sum(tz), np.sum(tc), np.sum(thc),np.sum(tcw))
    #print('before normalization tthin, tthick,tdisk,tbulge, all:', np.sum(tthin),  np.sum(tthick),  np.sum(tdisk), np.sum(tbulge), np.sum(tall))

    tthin =tthin *np.sum(wthin)/np.sum(tthin )
    tthick=tthick*np.sum(wthick)/np.sum(tthick)

    tdisk =tdisk *np.sum(wdisk)/np.sum(tdisk )
    tbulge=tbulge*np.sum(wbulge)/np.sum (tbulge)
    tall =tall *np.sum(wall)/np.sum(tall )

    #print( 'after normalization th  tw  tz tc thc tcware:', np.sum(th),  np.sum(tw), np.sum(tz), np.sum(tc), np.sum(thc), np.sum(tcw))
    #print( 'after normalization tthin, tthick,tdisk,tbulge:', np.sum(tthin),  np.sum(tthick),  np.sum(tdisk), np.sum(tbulge), np.sum(tall))

    #totalf = np.sum(th) + np.sum(tw) + np.sum(tz)+ np.sum(tc)

    ###if you want to check that the sum of the flux of the components is the same as the total, uncomment the two extra totalf and print them
    totalf = np.sum(tthin) + np.sum(tthick) + np.sum(tbulge)
    #totalf2 = np.sum(tdisk)+ np.sum(tbulge)
    #totalf3=np.sum (tall)
    #print('total fluxes with thin thick bulge', totalf,' and with disk bulge',totalf2, ' and tot', totalf3 )
    tthin =tthin /totalf
    tthick=tthick/totalf

    tdisk =tdisk /totalf
    tbulge=tbulge/totalf
    tall=tall/totalf
    #tcw=tcw/totalf
    #flux = th +  tw + tz + tc
    flux = tthin +  tthick + tbulge
    flux2= tdisk+tbulge
    #print("luminosity fractions f_thin, f_thick, f_disk, f_bulge, tot1, tot2 (disk+bulge), all")
    #print(np.sum(th),  np.sum(tw), np.sum(tz), np.sum(tc),np.sum(thc),np.sum(tcw))
    #print(np.sum(tthin),  np.sum(tthick),np.sum(tdisk), np.sum(tbulge), np.sum(flux), np.sum(flux2), np.sum(tall))





    ### SAVE DATA TO A FILE
   #if savedata:
   #    with open(figdir + 'SB_hmz.dat', 'w') as outfile:
   #        outfile.write('x/arcs,  y/arcs,  SB thin disk,  SB warm disk,  SB bulge, SB counter_rot' + '\n')
   #        for j in range(0, len(s)):
   #            outfile.write(("%10.4f" % xi[s[j]])  + ("%10.4f" % yi[s[j]]) +  ("%12.3e" % th[s[j]]) +
   #                          ("%12.3e" % tw[s[j]]) + ("%12.3e" % tz[s[j]]) + ("%12.3e" % tc[s[j]]) + '\n')






    hho_all = [ 'thin_d','warm_d', 'disk','bulge','all']


    ### EVALUATE VMAX and SMAX and SMIN

    vmax = np.nanmax([vel_thin,vel_thick, vel_disk, vel_bulge, vel_all])
    sig_t = np.array((sig_thin,sig_thick,sig_disk,sig_bulge, sig_all))
    #vmax = 79


    smax = np.nanmax(sig_t[sig_t > 0])
    smin = np.nanmin(sig_t[sig_t > 0])
    #smax =346
    #smin =161

    minf=min(-2.5 * np.log10(flux))
    maxf=max(-2.5 * np.log10(flux[flux !=0]))
    xi_t= (xi[s])
    yi_t=(yi[s])

    with open(wdir + str(gal)+'_kin_comps_test_s22_'+conversion+'.dat', 'w') as outfile:
            outfile.write('x/arcs  y/arcs  SB_thin_disk  SB_thick_disk  SB_disk  SB_bulge SB_whole vel_thin_disk  vel_thick_disk  vel_disk  vel_bugle vel_whole sig_thin_disk  sig_thick_disk sig_disk  sig_bulge sig_whole' + '\n')
            for j in range(0, len(s)):
                outfile.write(("%10.4f" % xi_t[j])  + ("%10.4f" % yi_t[j]) +
                              ("%12.3e" % tthin[s[j]]) + ("%12.3e" % tthick[s[j]]) + ("%12.3e" % tdisk[s[j]])+ ("%12.3e" % tbulge[s[j]]) + ("%12.3e" % tall[s[j]]) +
                              ("%12.3e" % vel_thin[grid[s[j]]]) + ("%12.3e" % vel_thick[grid[s[j]]])+ ("%12.3e" % vel_disk[grid[s[j]]]) + ("%12.3e" % vel_bulge[grid[s[j]]])+ ("%12.3e" % vel_all[grid[s[j]]]) +
                              ("%12.3e" % sig_thin[grid[s[j]]]) + ("%12.3e" % sig_thick[grid[s[j]]])+ ("%12.3e" % sig_disk[grid[s[j]]]) + ("%12.3e" % sig_bulge[grid[s[j]]])+ ("%12.3e" % sig_all[grid[s[j]]]) +'\n')

    # print ('VMAX, SMAX, SMIN are:', vmax, smax, smin)
    # print(np.max(th),np.min(th),np.max(tw),np.min(tw),np.max(tz),np.min(tz),np.max(tc),np.min(tc))


    ### PLOT THE RESULTS
    # Plot settings
    figfile=wdir+gal+'_comps_kin_test_s22_'+conversion+'.pdf'
    plt.figure(figsize=(12, 18))
    #plt.subplots_adjust(hspace=0.7, wspace=0.01, left=0.01, bottom=0.05, top=0.99, right=0.99)
    plt.subplots_adjust(hspace=0.4, wspace=0.02, left=0.01, bottom=0.05, top=0.99, right=0.99)




    ### PLOT THE COMPONENTS
    ## COLD
    ax1=plt.subplot(5, 3, 1)
    display_pixels(xi[s], yi[s], -2.5 * np.log10(tthin[s]) , pixelsize=dx,
                    colorbar=True, nticks=7, cmap='YlOrRd_r', label='-2.5 log10(flux)', vmin=minf, vmax=maxf)
    #ellipse = Ellipse((0, 0),width=Re * 2,height=semi_min * 2,facecolor='none', edgecolor = 'black')
   # ax1.add_patch(ellipse)
    ax2=plt.subplot(5, 3, 2)
    plt.title("THIN DISK COMPONENT")
    display_pixels(xi[s], yi[s], vel_thin[grid[s]], pixelsize=dx, colorbar=True, nticks=7, cmap='RdYlBu_r',
                   vmin=-1.0 * vmax, vmax=vmax, label='Velocity')

    #ellipse = Ellipse((0, 0),width=Re * 2,height=semi_min * 2,facecolor='none', edgecolor = 'black')
    #ax2.add_patch(ellipse)

    ax3=plt.subplot(5, 3, 3)
    display_pixels(xi[s], yi[s], sig_thin[grid[s]], pixelsize=dx, colorbar=True, nticks=7, cmap='YlOrRd',
                   vmin=smin, vmax=smax, label=r'$\mathbf{\sigma}$')
    #ellipse = Ellipse((0, 0),width=Re * 2,height=semi_min * 2,facecolor='none', edgecolor = 'black')
    #ax3.add_patch(ellipse)
    ## WARM
    plt.subplot(5, 3, 4)
    display_pixels(xi[s], yi[s], -2.5 * np.log10(tthick[s]) , pixelsize=dx,
                   colorbar=True, nticks=7, cmap='YlOrRd_r', label='-2.5 log10(flux)', vmin=minf, vmax=maxf)
    plt.subplot(5, 3, 5)
    plt.title("THICK DISK COMPONENT")
    display_pixels(xi[s], yi[s], vel_thick[grid[s]], pixelsize=dx, colorbar=True, nticks=7, cmap='RdYlBu_r',
                   vmin=-1.0 * vmax, vmax=vmax, label='Velocity')
    plt.subplot(5, 3, 6)
    display_pixels(xi[s], yi[s], sig_thick[grid[s]], pixelsize=dx, colorbar=True, nticks=7, cmap='YlOrRd',
                   vmin=smin, vmax=smax, label=r'$\mathbf{\sigma}$')

      #HOT+CR
    plt.subplot(5, 3, 7)
    display_pixels(xi[s], yi[s], -2.5 * np.log10(tdisk[s]) , pixelsize=dx,
                   colorbar=True, nticks=7, cmap='YlOrRd_r', label='-2.5 log10(flux)', vmin=minf, vmax=maxf)
    plt.subplot(5, 3, 8)
    #plt.title("HOT + C.R. COMPONENT")
    plt.title("DISK COMPONENT")
    display_pixels(xi[s], yi[s], vel_disk[grid[s]], pixelsize=dx, colorbar=True, nticks=7, cmap='RdYlBu_r',
                   vmin=-1.0 * vmax, vmax=vmax, label='Velocity')
    plt.subplot(5, 3, 9)
    display_pixels(xi[s], yi[s], sig_disk[grid[s]], pixelsize=dx, colorbar=True, nticks=7, cmap='YlOrRd',
                   vmin=smin, vmax=smax, label=r'$\mathbf{\sigma}$')


    plt.subplot(5, 3, 10)
    display_pixels(xi[s], yi[s], -2.5 * np.log10(tbulge[s]) , pixelsize=dx,
                   colorbar=True, nticks=7, cmap='YlOrRd_r', label='-2.5 log10(flux)', vmin=minf, vmax=maxf)
    plt.subplot(5, 3, 11)
    #plt.title("HOT + C.R. COMPONENT")
    plt.title("BULGE COMPONENT")
    display_pixels(xi[s], yi[s], vel_bulge[grid[s]], pixelsize=dx, colorbar=True, nticks=7, cmap='RdYlBu_r',
                   vmin=-1.0 * vmax, vmax=vmax, label='Velocity')
    plt.subplot(5, 3, 12)
    display_pixels(xi[s], yi[s], sig_bulge[grid[s]], pixelsize=dx, colorbar=True, nticks=7, cmap='YlOrRd',
                   vmin=smin, vmax=smax, label=r'$\mathbf{\sigma}$')
    #########all
    plt.subplot(5, 3, 13)
    display_pixels(xi[s], yi[s], -2.5 * np.log10(tall[s]) , pixelsize=dx,
                   colorbar=True, nticks=7, cmap='YlOrRd_r', label='-2.5 log10(flux)', vmin=minf, vmax=maxf)
    plt.subplot(5, 3, 14)
    #plt.title("HOT + C.R. COMPONENT")
    plt.title("WHOLE COMPONENT")
    display_pixels(xi[s], yi[s], vel_all[grid[s]], pixelsize=dx, colorbar=True, nticks=7, cmap='RdYlBu_r',
                   vmin=-1.0 * vmax, vmax=vmax, label='Velocity')
    plt.subplot(5, 3, 15)
    display_pixels(xi[s], yi[s], sig_all[grid[s]], pixelsize=dx, colorbar=True, nticks=7, cmap='YlOrRd',
                   vmin=smin, vmax=smax, label=r'$\mathbf{\sigma}$')


    plt.tight_layout()
    plt.savefig(figfile)
    plt.close()

    return 1



def run_dec(gal, w_dir1, Re):


    figdir=w_dir1 + gal + '/figure_nn/' # absolute path to keep figures
    best_fit= ascii.read(figdir+ gal+'_best_fit.dat', data_start=0,guess=False, delimiter=' ') #file with the bestfit values (these are also the name of the folders where the orbit libraries are)
    #print(best_fit[0][1])
    diri_o= str(best_fit[0][0])   # the subdirectory contains the best models
    ml_str=str(best_fit[0][1])                               # the best model
    print(diri_o, ml_str)

    rootname=''
    wdir=w_dir1+ str(gal) +'/'+diri_o +'/'
    bf_dir=wdir+ml_str
    aph_bf=bf_dir+'/nn_aphist.out'

    head1 = np.genfromtxt(aph_bf, max_rows=1)
    w = int(head1[0]) #bins
    n_apertures = int(head1[1]) #apertures (kinematic bins)

    ocut = [0.8, 0.25,-0.25] #cuts in lambda_z for the components (as in Santucci+22)



    Rmax_arcs=15#gal_infos['Rmax[arcsec]'][i]
    xrange = None

    scale_factor = scale_factor_find(bf_dir) #read the scale factor used to rescale the orbits

    #create_orbital_component_files(w_dir=w_dir1, object=gal, rootname=rootname, diri_o=diri_o, ml_str = ml_str,
    #                   ocut=ocut, Rmax_arcs=Rmax_arcs, xrange=xrange)


    losvd_histograms, intrinsic_masses, n_orbs, projected_masses=read_losvd_histograms(wdir, n_apertures,scale_factor)  #read the losvd of the orbits
    #print('losvd_histograms, intrinsic_masses, n_orbs, projected_masses')
    #print(losvd_histograms, intrinsic_masses, n_orbs, projected_masses)
    #comps_aphist(w_dir=w_dir1, object=object, aph_bf=aph_bf, vscale=scale_factor)

    #plot_comps(w_dir=w_dir1, object=object, rootname=None,savedata=True, xlim=15 , ylim=15)
    print('norbs',n_orbs)

    print('losvd shape', losvd_histograms.y.shape)
    create_orbital_component_files_giu(w_dir=w_dir1, object=gal, rootname=rootname, diri_o=diri_o, ml_str = ml_str,
                      ocut=ocut, Rmax_arcs=Rmax_arcs, xrange=xrange) #create the files with the orbits selected for each components

    return losvd_histograms, projected_masses




w_dir='/Users/z5178033/Downloads/' # the working directory where you have all the galaxy catalogue
w_dir='/Users/maindl/Documents/DriveM/drtim/ClientsNonGAMS/DYNAMITE/Giulia/decomposition/'

gal_infos=ascii.read(w_dir+'SAMI_gals_cat.dat')
gals=gal_infos['CATAID']




gal=41059
print('DOING GALAXY ', gal)

i = np.where(gal_infos['CATAID'] == gal)[0][0]
Re=gal_infos['Re[arcs]'][i] #
w_dir1='/Users/z5178033/Downloads/test_crcut/' #folder where you have the galaxies folders with the data, if different from where you have the catalogue
w_dir1='/Users/maindl/Documents/DriveM/drtim/ClientsNonGAMS/DYNAMITE/Giulia/decomposition/'
losvd_histograms, proj_mass=run_dec(str(gal), w_dir1, Re) #read the orbits and create the velocity histrogram
print('Orbits read')
comps=['disk', 'thin_d', 'warm_d', 'bulge', 'all']




# for gal in gals:
#     if gal==41059:
#         print('DOING GALAXY ', gal)
#
#         i = np.where(gal_infos['CATAID'] == gal)[0][0]
#         Re=gal_infos['Re[arcs]'][i] #
#         w_dir1='/Users/z5178033/Downloads/test_crcut/' #folder where you have the galaxies folders with the data, if different from where you have the catalogue
#         losvd_histograms, proj_mass=run_dec(str(gal), w_dir1, Re) #read the orbits and create the velocity histrogram
#         print('Orbits read')
#         comps=['disk', 'thin_d', 'warm_d', 'bulge', 'all']
#         conversion='gh_fit'#other options are 'losvd_vsig', 'fortran', 'moments'
#         comps_aphist(w_dir1, losvd_histograms, comps, str(gal), conversion)#select the components and calculate the kinematics for each (this is done with the selection used in Santucci+22)
#         print('Compsonents done')
#         plot_comps_giu(w_dir=w_dir1, gal=str(gal), rootname=None,savedata=True, xlim=15 , ylim=15, Re=Re) #plot the kinematics
#         print('Plots done')
#         print('**************************')
#



conversion='gh_fit'#other options are 'losvd_vsig', 'fortran', 'moments'
comps_aphist(w_dir1, losvd_histograms, comps, str(gal), conversion)#select the components and calculate the kinematics for each (this is done with the selection used in Santucci+22)
print('Compsonents done')
plot_comps_giu(w_dir=w_dir1, gal=str(gal), rootname=None,savedata=True, xlim=15 , ylim=15, Re=Re, conversion=conversion) #plot the kinematics
print('Plots done')
print('**************************')




conversion='losvd_vsig'#other options are 'losvd_vsig', 'fortran', 'moments'
comps_aphist(w_dir1, losvd_histograms, comps, str(gal), conversion)#select the components and calculate the kinematics for each (this is done with the selection used in Santucci+22)
print('Compsonents done')
plot_comps_giu(w_dir=w_dir1, gal=str(gal), rootname=None,savedata=True, xlim=15 , ylim=15, Re=Re, conversion=conversion) #plot the kinematics
print('Plots done')
print('**************************')




conversion='fortran'#other options are 'losvd_vsig', 'fortran', 'moments'
comps_aphist(w_dir1, losvd_histograms, comps, str(gal), conversion)#select the components and calculate the kinematics for each (this is done with the selection used in Santucci+22)
print('Compsonents done')
plot_comps_giu(w_dir=w_dir1, gal=str(gal), rootname=None,savedata=True, xlim=15 , ylim=15, Re=Re, conversion=conversion) #plot the kinematics
print('Plots done')
print('**************************')




conversion='moments'#other options are 'losvd_vsig', 'fortran', 'moments'
comps_aphist(w_dir1, losvd_histograms, comps, str(gal), conversion)#select the components and calculate the kinematics for each (this is done with the selection used in Santucci+22)
print('Compsonents done')
plot_comps_giu(w_dir=w_dir1, gal=str(gal), rootname=None,savedata=True, xlim=15 , ylim=15, Re=Re, conversion=conversion) #plot the kinematics
print('Plots done')
print('**************************')

