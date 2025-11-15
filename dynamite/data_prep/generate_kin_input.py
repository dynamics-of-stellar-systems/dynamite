#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2021-2024

@author: Sabine Thater, Julia Lamprecht, Thomas I. Maindl
"""
import numpy as np
from astropy import table
from astropy.io import fits
import matplotlib.pyplot as plt
import plotbin.display_pixels as dp
import pafit.fit_kinematic_pa as pa
import cmasher as cmr


def read_kinematics_user(file):
    # This function can be filled by the user
    return None, None, None, None


def read_kinematics_nifs(f_names, n_gh='all', idl=False):
    """
    Read NIFS Gauss Hermite kinematics data and Voronoi binning

    Parameters
    ----------
    f_names : tuple or list of str
        The first item is the filename of the file holding kinematics data,
        the second is the filename of the Voronoi binning data.
        Format of f_names[0]: bin_number, v[km/s], dv[km/s],
        sigma[km/s], dsigma[km/s], #GH_moments (ignored), h3, dh3, h4, dh4,...
        Format of f_names[1]: x[arcsec], y[arcsec], bin_number
        Both: lines starting with # are ignored
    n_gh : int or str, optional
        If int: number of GH moments to read, will be capped at the number of
        moments provided in the data file. Must be >= 3. If 'all': all GH
        moments in the data file will be read. The default is 'all'.
    idl : bool, optional
        If True, bin_number is replaced by bin_number - 1 after reading the
        binning file.
        This is necessary if the binning file was created using IDL (in that
        case the bin_number in the file starts with 1). The default is False.

    Returns
    -------
    bin_num : 1d numpy array (int)
        Spatial bin ids.
    x : 1d numpy array (float)
        x coordinates mapped to bin id.
    y : 1d numpy array (float)
        y coordinates mapped to bin id..
    data_table : astropy table
        Kinematics. Columns: v, dv, sigma, dsigma, h3, dh3...h[n_gh], dh[n_gh].

    """
    # Load data from file
    data = np.loadtxt(f_names[0])  # kinematics data
    # Extract relevant columns from the data
    data = np.delete(data, [0,5], axis=1)  # ignore binNum, #GH moments columns
    highest_gh_data = data.shape[1] // 2  # highest GH moment in the data
    if n_gh == 'all':
        n_gh = highest_gh_data
    else:
        if n_gh > highest_gh_data:
            n_gh = highest_gh_data
            print(f'WARNING: {n_gh=} too big, reset to number of GH moments '
                  f'in data: {highest_gh_data}')
    # delete unused GH moments from data
    gh_to_delete = range(n_gh * 2, data.shape[1])
    data = np.delete(data, gh_to_delete, axis=1)
    # Note: Assuming specific column indices for different kinematic parameters
    # Adjust these indices based on the structure of the data files
    data_table = table.Table(data,names=['v','dv','sigma','dsigma'] +
                                        [f'{s}{i}' for i in range(3,n_gh + 1)
                                                   for s in ['h','dh']])
    # Subtract the median from the odd moments (including v)
    data_table['v'] -= np.median(data_table['v'])
    for i_gh in range(3, n_gh + 1, 2):
        data_table[f'h{i_gh}'] -= np.median(data_table[f'h{i_gh}'])

    # voronoi = np.loadtxt(f_names[1])  # Voronoi 2d binning
    # Extract bin numbers and spatial coordinates
    x, y, bin_num = np.loadtxt(f_names[1], unpack=True)  # Voronoi 2d binning
    bin_num = bin_num.astype(int)
    if idl:  # IDL starts at 1, python at 0
        bin_num -= 1

    return bin_num, x, y, data_table


def read_atlas3d(file):
    #adapted from http://www-astro.physics.ox.ac.uk/atlas3d/
    #read original cube

    hdu = fits.open(file[0])
    spectrum = hdu[0].data
    table = hdu[2].data
    #hdr = hdu[2].header #unused

    x = table["A"] # Coordinates of the original spaxels in arcsec
    y = table["D"]
    flux = np.mean(spectrum, 1)  #surface brightness

    hdu = fits.open(file[1])
    table = hdu[1].data
    #kin_hdr=hdu[1].header #unused

    xgen = table['XS'] # Voronoi generators
    ygen = table['YS']
    vel = table['VPXF']
    sig=table['SPXF']
    h3=table['H3PXF']
    h4=table['H4PXF']
    dvel = table['EVPXF']
    dsig=table['ESPXF']
    dh3=table['EH3PXF']
    dh4=table['EH4PXF']

    #perform Voronoi tesselation starting from the nodes values, adapted from axisymm Schwarzschild code
    npixels=len(x)
    binNum=np.zeros((npixels), dtype=int)
    for j, (xj, yj) in enumerate(zip(x, y)):
        binNum[j] = np.argmin(((xj - xgen)**2 + (yj - ygen)**2))

    return binNum,x,y,flux,vel,sig,h3,h4,dvel,dsig,dh3,dh4,xgen,ygen
    #can be switched on for a check. They should show the same. If not, something went wrong with the V. tesselation
    #plt.figure()
    #im1=dbg.display_bins_generators(xgen,ygen,vel,x,y,label='Velocity (km/s)')
    #plt.figure()
    #im1=dp.display_pixels(x,y,vel[binNum],label='Velocity (km/s)')


def read_califa(file):
    hdulist = fits.open(file)
    kin_tab = hdulist[1].data
    #kin_hdr = hdulist[1].header #unused

    s = kin_tab.BIN_ID > 0

    binNum = (kin_tab[s].BIN_ID).astype(int)

    xp = kin_tab[s].X
    yp = kin_tab[s].Y
    flux=kin_tab[s].FLUX

    #check where each bin appears the first time
    ubins, indices = np.unique(binNum, return_index=True)

    vel = kin_tab[s].VPXF[indices]
    sig = kin_tab[s].SPXF[indices]
    h3 = kin_tab[s].H3PXF[indices]
    h4 = kin_tab[s].H4PXF[indices]
    dvel = kin_tab[s].DVPXF[indices]
    dsig = kin_tab[s].DSPXF[indices]
    dh3 = kin_tab[s].DH3PXF[indices]
    dh4 = kin_tab[s].DH4PXF[indices]

    # adopted from schwpy routine
    # correct on minium error on higher moment values from the ppxf paper
    minvelerr = 1.0
    dvel[dvel < minvelerr] = minvelerr
    dsig[dsig < minvelerr] = minvelerr
    mingherr = 0.005
    dh3[dh3 < mingherr] = mingherr
    dh4[dh4 < mingherr] = mingherr

    # adopted from schwpy routine
    # Check that the dispersion is positive
    if np.min(sig) < 0:
        print('Some dispersion values are negative! Correcting!')
        mask = sig < 0
        sig[mask] = np.median(sig[sig >= 0])
        dsig[mask] = 1e5

    nbin = int(np.max(kin_tab.BIN_ID))
    xbin = np.zeros(nbin)
    ybin = np.zeros(nbin)
    for i in range(len(xbin)):
        si = np.where(ubins == i+1)[0]
        xbin[i] = kin_tab[si].X
        ybin[i] = kin_tab[si].Y

    return binNum-1,xp,yp,flux,vel,sig,h3,h4,dvel,dsig,dh3,dh4,xbin,ybin


def create_aperture_file(dir,expr,minx,maxx,miny,maxy,angle_deg,nx,ny):
    # The angle saved here is measured counter clock-wise
    # from the galaxy major axis to the X-axis of the input data
    aperture_file = open(dir+'aperture'+expr+'.dat', 'w')
    # The following line is a comment and optional.
    aperture_file.write('#counter_rotation_boxed_aperturefile_version_2 \n')
    # Now, write the actual data...
    aperture_file.write('\t{0:<.6f}\t{1:<.6f} \n'.format(minx, miny))
    aperture_file.write('\t{0:<.6f}\t{1:<.6f} \n'.format(maxx-minx, maxy-miny))
    aperture_file.write('\t{0:<.6f} \n'.format(angle_deg))
    aperture_file.write('\t{0}\t{1} \n'.format(int(nx), int(ny)))
    aperture_file.close()


def create_bins_file(dir,expr,grid):
    #adapted from schwpy. Not very beautiful

    s = np.shape(grid)
    #the grid is expressed with 10 numbers in each line. This is prepared here
    flattened = grid.T.flatten()
    num_of_lines = int(len(flattened)/10)
    bins_file = open(dir+'bins'+expr+'.dat', 'w')
    # The following line is a comment and optional.
    bins_file.write('#Counterrotation_binning_version_1\n')
    # Now, write the actual data...
    bins_file.write('{0}\n'.format(int(s[0]*s[1])))

    for line in range(num_of_lines):
        string = '\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(int(flattened[line*10]),int(flattened[line*10+1]),
                                                                               int(flattened[line*10+2]),int(flattened[line*10+3]),
                                                                               int(flattened[line*10+4]),int(flattened[line*10+5]),
                                                                               int(flattened[line*10+6]),int(flattened[line*10+7]),
                                                                               int(flattened[line*10+8]),int(flattened[line*10+9]))
        bins_file.write(string)
    last_line = ''
    # number of items in last line:
    num_last_items = len(flattened) - num_of_lines*10
    for i in range(num_last_items):
        last_line += '\t{0}'.format(int(flattened[i+num_of_lines*10]))
    last_line += '\n'
    bins_file.write(last_line)
    bins_file.close()


def kin_file(directory,expr,data):
    data.write(directory + 'gauss_hermite_kins' + expr + '.ecsv',
               format='ascii.ecsv',
               overwrite=True)


#### main routine ######
def create_kin_input(galaxy, file, dyn_model_dir, expr='', angle_deg=0,
                     ngh='all',
                     pointsym=0, bisym=0, xoffset=0, yoffset=0, voffset=0,
                     fit_PA=False, kin_input=0, files=True, plot=True):
    """
    Create the inputs from the Voronoi binned kinematics files
    Returns bins.dat, aperture.dat and kin_data.dat and plots the kinematics

    Note: sym keywords will be added later!
    """
    print('Galaxy: {0}'.format(galaxy))
    print(file)

    #read in binned kinematics, this needs to be changed by the user
    if kin_input == 'CALIFA':
        if ngh == 'all':
            ngh = 4
        binNum,xp,yp,flux,vel,sig,h3,h4,dvel,dsig,dh3,dh4,xbin,ybin=read_califa(file)

    elif kin_input == 'ATLAS3D':
        if ngh == 'all':
            ngh = 4
        binNum,xp,yp,flux,vel,sig,h3,h4,dvel,dsig,dh3,dh4,xbin,ybin=read_atlas3d(file)

    elif kin_input == 'USER':
        binNum, xp, yp, data = read_kinematics_user(file)

    elif kin_input == 'NIFS':
        binNum, xp, yp, data = read_kinematics_nifs(file, n_gh=ngh, idl=True)

    if kin_input in ['CALIFA', 'ATLAS3D']:
        data = table.Table()
        data['v'] = vel - np.median(vel)  # not done in califa and atlas3d read
        data['dv'] = dvel
        data['sigma'] = sig
        data['dsigma'] = dsig
        data['h3'] = h3
        data['dh3'] = dh3
        data['h4'] = h4
        data['dh4'] = dh4
        if ngh == 6:
            data['h5'] = np.full_like(h4, 0)
            data['dh5'] = np.full_like(h4, 0.3)
            data['h6'] = np.full_like(h4, 0)
            data['dh6'] = np.full_like(h4, 0.3)
    n_gh = len(data.colnames) // 2
    # add first column holding a bin index starting at 1
    data.add_column(np.arange(1, len(data) + 1), name='vbin_id', index=0)

    xp=xp+xoffset
    yp=yp+yoffset
    # xbin=xbin+xoffset  # not used
    # ybin=ybin+yoffset  # not used

    map1 = cmr.get_sub_cmap('twilight_shifted', 0.05, 0.6)
    map2 = cmr.get_sub_cmap('twilight_shifted', 0.05, 0.95)

    if fit_PA:
        angle_deg,_,vel_syst = pa.fit_kinematic_pa(xp,
                                                   yp,
                                                   data['v'][binNum],
                                                   cmap=map2)
        plt.savefig(dyn_model_dir + 'pafit' + expr + '.pdf')
        data['v'] -= vel_syst
    data['v'] += voffset  # add velocity offset if required

    # Determination of the pixel size
    npixels = len(xp)
    dx = 1e30
    for j in range(npixels-1):
        dxj = np.min((xp[j]-xp[j+1])**2 + (yp[j]-yp[j+1])**2)
        if dxj < dx:
            dx = dxj
    dx = np.sqrt(dx)

    print('The pixel size is {0}'.format(dx))
    # Adjust position angle
    print('PA: {0}'.format(angle_deg))
    angle_deg = 90 - angle_deg
    print('Total bins: {0}'.format(len(data)))

    maxx = np.max(xp) + dx/2.0
    maxy = np.max(yp) + dx/2.0
    minx = np.min(xp) - dx/2.0
    miny = np.min(yp) - dx/2.0
    nx = np.round((maxx-minx)/dx)
    ny = np.round((maxy-miny)/dx)
    print("Pixel grid dimension is ", nx, ny)

    grid = np.zeros((int(nx), int(ny)))
    k = ((xp-minx)/dx).astype(int)
    j = ((yp-miny)/dx).astype(int)
    grid[k, j] = binNum+1

    if plot is True:
        # Take only 98 / 2 for plotting to get rid of outliers
        vmax = np.percentile(np.abs(data['v']), 98)
        smax = np.percentile(data['sigma'], 98)
        smin = np.percentile(data['sigma'], 2)
        print('Vels plot: {0}, {1}, {2}'.format(vmax, smin, smax))

        plot_rows = (n_gh - 1) // 6 + 1
        if n_gh <= 4:
            plot_cols = 4
            figsize = (18, 4)
        else:
            plot_cols = 6
            figsize = (25, 5 * plot_rows)
        fig, ax = plt.subplots(plot_rows, plot_cols, figsize=figsize)
        if n_gh % plot_cols > 0:  # remove axes from empty subplots
            for col in range(n_gh % plot_cols, plot_cols):
                if plot_rows == 1:
                    ax[col].set_axis_off()
                else:
                    ax[plot_rows - 1, col].set_axis_off()

        plt.subplot(plot_rows,plot_cols,1)
        plt.title('Velocity [km/s]')
        plt.xlabel('arcsec', fontsize=10)
        plt.ylabel('arcsec', fontsize=10)
        dp.display_pixels(xp, yp, data['v'][binNum],
                          angle=angle_deg,
                          vmin=-vmax, vmax=vmax,
                          cmap=map2)

        plt.subplot(plot_rows,plot_cols,2)
        plt.title('Velocity dispersion [km/s]')
        plt.xlabel('arcsec', fontsize=10)
        dp.display_pixels(xp, yp, data['sigma'][binNum],
                          angle=angle_deg,
                          vmin=smin, vmax=smax,
                          cmap=map1)

        for h_mom in range(3, n_gh + 1):
            plt.subplot(plot_rows,plot_cols,h_mom)
            plt.title(f'$h_{{{h_mom}}}$ moment')
            plt.xlabel('arcsec', fontsize=10)
            if (h_mom - 1) % plot_cols == 0:
                plt.ylabel('arcsec', fontsize=10)
            dp.display_pixels(xp, yp, data[f'h{h_mom}'][binNum],
                              angle=angle_deg,
                              vmin=-0.15, vmax=0.15,
                              cmap=map2)

        fig.savefig(dyn_model_dir + 'kinmaps' + expr + '.pdf')

    # round v, dv, sigma, dsigma, and the GH moments (but not their errors)
    data.round(decimals={c:4 for c in ['v','dv','sigma','dsigma'] +
                                      [f'h{i}' for i in range(3,n_gh + 1)]})

    if files is True:
        create_aperture_file(dyn_model_dir,
                             expr,
                             minx, maxx, miny, maxy,
                             angle_deg,
                             nx, ny)
        create_bins_file(dyn_model_dir, expr, grid)
        kin_file(dyn_model_dir, expr, data)
