#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 16:10:16 2025

@author: ggslc
"""
import numpy as np
from downscale import interp_to_local, interp_to_surface
from transformation import local_to_global_map, cell_id


def glint_conservation_adjust(smb_ism_nc, mask_ism, grid_ism,
                              area_atm, grid_atm, valid_atm,
                              down_transform):

    # conservation enforcement - see glint_downscaling_gcm
    # roughly - adjust smb so that the total smb_ism within an
    # atmosphere grid box equals the smb_atm integrated vertically.
    # this resembles the glint code. blocky.
    # there is a NASTY! comment in the glint code
    smb_ism = np.zeros(smb_ism_nc.shape) + smb_ism_nc

    ism_to_atm_map = local_to_global_map(down_transform, grid_atm, grid_ism)

    dx = grid_ism.axes[0][1] - grid_ism.axes[0][0]  # \todo - move to Uniform2DGrid

    smb_adjust = np.zeros(grid_atm.shape)
    for I in grid_atm.axes_index[0]:
        for J in grid_atm.axes_index[1]:
            if np.any(valid_atm[:, J, I]):

                cc = cell_id(I, J, *grid_atm.shape)
                
                # vertical sum
                atm_col_total = np.sum(smb_atm[:, J, I]*area_atm[:, J, I])
    
                # horizontal sum
                ism_total = dx**2 * np.sum(np.where(ism_to_atm_map == cc, smb_ism, 0.0))
                ism_area =  dx**2 * np.sum(np.where(ism_to_atm_map == cc, 1.0, 0.0))

                if ism_area > 0:
                    smb_adjust[I, J] = (atm_col_total - ism_total)/ism_area
                    smb_ism = np.where(ism_to_atm_map == cc,
                                       smb_ism + smb_adjust[I, J],
                                       smb_ism)

    return smb_ism * mask_ism


def atm_to_ism(smb_atm, sfct_atm, topo_atm, area_atm,  grid_atm,
               topo_ism, frac_ism,  grid_ism,
               up_transform,  down_transform):

    atm_coords = down_transform.local_xy((grid_atm.coords))

    valid_atm = area_atm.data > 0.0
    smb_xyz = interp_to_local(smb_atm.data, atm_coords, grid_ism.coords,
                              order=1, valid=valid_atm & ~smb_atm.mask)

    sfct_xyz = interp_to_local(sfct_atm.data, atm_coords, grid_ism.coords,
                               order=1, valid=~sfct_atm.mask)

    # this is probably uniform in x,y , but ...
    topo_xyz = interp_to_local(topo_atm, atm_coords, grid_ism.coords, order=1)

    mask_ism = np.where(frac_ism > 0.01, 1, 0)

    sfct_ism, smb_ism = \
        (mask_ism*interp_to_surface(f_xyz, topo_xyz, topo_ism)
         for f_xyz in (sfct_xyz, smb_xyz))

    if not type(area_atm) is None:
        smb_ism = glint_conservation_adjust(
            smb_ism, mask_ism, grid_ism,
            area_atm, grid_atm, valid_atm,
            down_transform)

    return smb_ism, sfct_ism


def ism_to_atm():

    icefrac_atm = None

    return icefrac_atm


if __name__ == "__main__":

    #test program - move to test_coupling.py
    from netCDF4 import Dataset
    from transformation import up_down_pair, Uniform2DGrid, Box, PROJ_ARCTIC_4326
    import matplotlib.pyplot as plt

    #transformation pair
    up_tr, down_tr = up_down_pair(PROJ_ARCTIC_4326)

    #example ice sheet model output (input to the coupler). bisicles CF
    nc_ism_in = Dataset('bike_cf_gris4326.nc','r')
    topo_ism = nc_ism_in['orog'][:,:]
    frac_ism = nc_ism_in['sftgif'][:,:]
    grid_ism = Uniform2DGrid(nc_ism_in['x'][:].data, nc_ism_in['y'][:].data)


    #example atmosphere model output (input to the coupler).
    nc_atm_in = Dataset('atmos_cx209c_P1Y_20000101-20010101_icecouple.nc','r')
    lon_ism, lat_ism = up_tr.global_XY(grid_ism.coords)
    #plt.pcolormesh(lon_ism, lat_ism,topo_ism,vmin=0,vmax=2000, cmap='jet')


    ilo, ihi = 148, 188
    #ilo, ihi = 0, -1
    jlo, jhi = 116, 140

    #ilo, ihi = 0, 192
    #jlo, jhi = 0, 144
    grid_atm = Uniform2DGrid(nc_atm_in['longitude'][ilo:ihi].data,
                                    nc_atm_in['latitude'][jlo:jhi].data)

    prep_atm = lambda arr: np.transpose(arr[:,jlo:jhi,ilo:ihi], axes=[0,1,2])

    area_atm = prep_atm(nc_atm_in['tile_surface_area'])
    smb_atm = prep_atm(nc_atm_in['ice_smb']) * 910 * 12 # values suggest kg/month
    sfct_atm = prep_atm(nc_atm_in['ice_stemp'])

    lon, lat = grid_atm.coords
    nlat, nlon = lat.shape
    nec =  10
    topo_atm = np.zeros(smb_atm.shape)
    #making this up too
    topo_mid = np.array([31.86, 297.01, 551.87, 846.16, 1151.70,
                1457.07, 1808.83, 2257.02,  2737.89, 3099.39])
    for ec in range(0,nec):
        topo_atm[ec,:,:] = topo_mid[ec]


    plt.pcolormesh(*grid_atm.axes, np.mean(sfct_atm, axis=0),vmin=-50, cmap='bwr_r')
    plt.colorbar()
    plt.show()

    #raise



#%% downscaling & upscaling
    smb_ism, sfct_ism = atm_to_ism(smb_atm, sfct_atm, topo_atm, area_atm,
                                   grid_atm, topo_ism, frac_ism,
                                   grid_ism, up_tr, down_tr)
# %% plot
    if True:

        def cf(ax, z, zl, label=True):
            cs = ax.contour(*km(*grid_ism.axes), z, zl, colors=['k'] ,linewidths=0.5, linestyles='-')
            if label:
                ax.clabel(cs, cs.levels, fontsize=8)

        fig, axs = plt.subplots(1,2)
        [ax.set_aspect('equal') for ax in axs.flat]
        km = lambda x,y : (x*1.0e-3,y*1e-3)
        #a = 0.1
        ax = axs.flat[0]
        pc = ax.pcolormesh(*km(*grid_ism.axes),smb_ism,cmap='bwr_r',vmin=-1,vmax=1)
        fig.colorbar(pc,ax=ax,shrink=0.5)
        cf(ax, lon_ism, 360+np.arange(-60,-15,5))
        cf(ax, lat_ism, np.arange(60,90,5))
        cf(ax, frac_ism, [0.05], label=False)

        ax = axs.flat[1]
        pc = ax.pcolormesh(*km(*grid_ism.axes),sfct_ism,cmap='plasma',vmin=-50,vmax=0)
        fig.colorbar(pc,ax=ax,shrink=0.5)
        cf(ax, lon_ism, 360+np.arange(-60,-15,5))
        cf(ax, lat_ism, np.arange(60,90,5))
        cf(ax, frac_ism, [0.05], label=False)

