#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 16 21:39:16 2021

@author: sabine
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import plotbin.display_pixels as dp
import pafit.fit_kinematic_pa as pa
from astropy import table
import cmasher as cmr

def read_kinematics_user(file):
    #This function can be filled by the user
    return None

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
    # The angled saved here is measured counter clock-wise
    # from the galaxy major axis to the X-axis of the input data
    aperture_file = open(dir+'aperture'+expr+'.dat', 'w')
    aperture_file.write('#counter_rotation_boxed_aperturefile_version_2 \n')
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
    bins_file.write('#Counterrotaton_binning_version_1\n')
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

def kin_file(dir,expr,data):
    data.write(dir+'gauss_hermite_kins'+expr +'.ecsv', format='ascii.ecsv', overwrite=True)

#### main routine ######
def create_kin_input(object, file, dyn_model_dir, expr='', angle_deg=0, ngh=4,
                     pointsym=0, bisym=0, xoffset=0, yoffset=0, voffset=0, fit_PA=False, kin_input=0, files=True,plot=True):

    """
    Create the inputs from the Voronoi binned kinematics files
    Returns bins.dat, aperture.dat and kin_data.dat and plots the kinematics

    Note: sym keywords will be added later!
    """
    print('Galaxy: {0}'.format(object))
    print(file)

    #read in binned kinematics, this needs to be changed by the user
    if kin_input == 'CALIFA':
        binNum,xp,yp,flux,vel,sig,h3,h4,dvel,dsig,dh3,dh4,xbin,ybin=read_califa(file)

    elif kin_input == 'ATLAS3D':
        binNum,xp,yp,flux,vel,sig,h3,h4,dvel,dsig,dh3,dh4,xbin,ybin=read_atlas3d(file)

    elif kin_input == 'USER':
        binNum,xp,yp,flux,vel,sig,h3,h4,dvel,dsig,dh3,dh4,xbin,ybin=read_kinematics_user(file)

    nbins = len(vel)

    xp=xp+xoffset
    yp=yp+yoffset
    xbin=xbin+xoffset
    ybin=ybin+yoffset

    map1 = cmr.get_sub_cmap('twilight_shifted', 0.05, 0.6)
    map2 = cmr.get_sub_cmap('twilight_shifted', 0.05, 0.95)

    vel = vel - np.median(vel)
    if fit_PA:
        angle_deg,_,vel_syst = pa.fit_kinematic_pa(xp,yp,vel[binNum],cmap=map2) 
        plt.savefig(dyn_model_dir+'pafit.pdf')
        vel=vel - vel_syst
    vel = vel + voffset

    # Determination of the pixel size
    npixels = len(xp)
    dx = 1e30
    for j in range(npixels-1):
        dxj = np.min((xp[j]-xp[j+1])**2 + (yp[j]-yp[j+1])**2)
        if dxj < dx:
            dx = dxj
    dx = np.sqrt(dx)

    print('The pixel size is {0}'.format(dx))
    print('PA: {0}'.format(angle_deg))
    angle_deg = 90 - angle_deg
    print('Total bins: {0}'.format(nbins))

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
        vmax = np.percentile(np.abs(vel), 98)
        smax = np.percentile(sig, 98)
        smin = np.percentile(sig, 2)
        print('Vels plot: {0}, {1}, {2}'.format(vmax, smin, smax))

        fig, axs = plt.subplots(1, 4, figsize=(18,4))
        plt.subplot(1,4,1)
        plt.title('Velocity [km/s]')
        dp.display_pixels(xp,yp,vel[binNum],angle=angle_deg,vmin=-vmax,vmax=vmax,cmap=map2)

        plt.subplot(1,4,2)
        plt.title('Velocity dispersion [km/s]')
        dp.display_pixels(xp,yp,sig[binNum],angle=angle_deg,vmin=smin,vmax=smax,cmap=map1)

        plt.subplot(1,4,3)
        plt.title(r'$h_{3}$ moment')
        dp.display_pixels(xp,yp,h3[binNum],angle=angle_deg,vmin=-0.15,vmax=0.15,cmap=map2)

        plt.subplot(1,4,4)
        plt.title(r'$h_{4}$ moment')
        dp.display_pixels(xp,yp,h4[binNum],angle=angle_deg,vmin=-0.15,vmax=0.15,cmap=map2)

        fig.subplots_adjust(left=0.04, wspace=0.3, hspace=0.01, right=0.97)
        fig.savefig(dyn_model_dir+'kinmaps.pdf')

    data = table.Table()
    data['vbin_id']=np.arange(1,nbins+1)
    data['v']=np.round(vel,decimals=4)
    data['dv']=np.round(dvel,decimals=4)
    data['sigma']=np.round(sig,decimals=4)
    data['dsigma']=np.round(dsig,decimals=4)
    data['h3']=np.round(h3,decimals=4)
    data['dh3']=dh3
    data['h4']=np.round(h4,decimals=4)
    data['dh4']=dh4

    if ngh==6:
        if (kin_input=='CALIFA') or ((kin_input=='ATLAS3D')):
            h5=np.full_like(h4, 0)
            h6=np.full_like(h4, 0)
            dh5=np.full_like(h4, 0.3)
            dh6=np.full_like(h4, 0.3)
        data['h5']=np.round(h5,decimals=4)
        data['dh5']=dh5
        data['h6']=np.round(h6,decimals=4)
        data['dh6']=dh6

    if files is True:
        create_aperture_file(dyn_model_dir,expr,minx,maxx,miny,maxy,angle_deg,nx,ny)
        create_bins_file(dyn_model_dir,expr,grid)
        kin_file(dyn_model_dir,expr,data)
