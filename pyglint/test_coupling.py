#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 09:48:10 2025

@author: ggslc
"""

import numpy as np
from coupling import atm_to_ism, ism_to_atm
from transformation import up_down_pair, Uniform2DGrid, PROJ_ARCTIC_4326
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from time import time

#%%
if __name__ == "__main__":


    def read_ism_nc(nc_file):
        nc_bike_in = Dataset(nc_file,'r')
        topo = nc_bike_in['orog'][:,:]
        frac = nc_bike_in['sftgif'][:,:]
        grid = Uniform2DGrid(nc_bike_in['x'][:].data, 
                             nc_bike_in['y'][:].data)
        mask = np.where(frac > 1e-4, True, False)
        
        return grid, topo, frac, mask

    def read_atm_nc(nc_file):
        
        nc_um_in = Dataset('atmos_cx209c_P1Y_20000101-20010101_icecouple.nc','r')
        ilo, ihi = 148, 188
        jlo, jhi = 119, 140
        ilo, ihi = 0, -1
        jlo, jhi = 0, -1
        grid = Uniform2DGrid(nc_um_in['longitude'][ilo:ihi].data,
                                nc_um_in['latitude'][jlo:jhi].data)

        def prep_um(arr):
            return arr[:,jlo:jhi,ilo:ihi]
        
        area = prep_um(nc_um_in['tile_surface_area'])
        smb = prep_um(nc_um_in['ice_smb']) * 910 * 12 # values suggest kg/month
        sfct = prep_um(nc_um_in['ice_stemp'])
        z_id = nc_um_in['tile_id']

        lon, lat = grid.coords
        nlat, nlon = lat.shape
        nec = z_id.shape[0]
        
        
        #making this up too
        topo_mid = np.array([31.86, 297.01, 551.87, 846.16, 1151.70,
                    1457.07, 1808.83, 2257.02,  2737.89, 3099.39])
        topo_max = np.zeros(nec+1) 
        topo_max[1:nec] = 0.5*(topo_mid[0:nec-1] + topo_mid[1:nec])
        topo_max[nec] = 5000.0
        topo = np.zeros(smb.shape)
        for ec in range(0,nec):
            topo[ec,:,:] = topo_mid[ec]  
            
        return grid, area, topo, topo_max, smb, sfct 
            

#%%
    #transformation pair
    up_tr, down_tr = up_down_pair(PROJ_ARCTIC_4326)

    #example ice sheet model output (input to the coupler). bisicles CF
    grid_bike, topo_bike, frac_bike, mask_bike = \
        read_ism_nc('bike_cf_gris4326.nc')
    

    #example atmosphere model output (input to the coupler).
    grid_um, area_um, topo_um, topo_max, smb_um, sfct_um = \
        read_atm_nc('atmos_cx209c_P1Y_20000101-20010101_icecouple.nc')
    

#%% downscaling & upscaling

    t = time()
    
    smb_bike, sfct_bike = atm_to_ism(smb_um, sfct_um, topo_um, area_um, grid_um, 
                                     topo_bike, frac_bike, mask_bike, grid_bike, 
                                     up_tr, down_tr)

    print(time() - t)
    t = time()
    
    um_shape = smb_um.shape
    frac_um = np.full(um_shape, 0.0)
    frac_um, topo_um = ism_to_atm(topo_max, area_um, grid_um,  
                         frac_um, topo_um,
                         topo_bike, frac_bike, mask_bike,  grid_bike,
                         up_tr, down_tr)
    
    print(time() - t)
# %% plot
    plot = True
    if plot:

        def cf(ax, z, zl, label=True):
            cs = ax.contour(*km(*grid_bike.axes), z, zl, colors=['k'] ,linewidths=0.5, linestyles='-')
            if label:
                ax.clabel(cs, cs.levels, fontsize=8)

        def km(x,y):
            return (x*1.0e-3,y*1e-3)

        lon_bike, lat_bike = up_tr(*grid_bike.coords)

        fig, axs = plt.subplots(2,2, figsize=(8,12))
        #__ = [ax.set_aspect('equal') for ax in axs.flat]
        #a = 0.1
        ax = axs.flat[0]
        pc = ax.pcolormesh(*grid_um.axes, np.mean(sfct_um, axis=0),vmin=-50, cmap='bwr_r')
        fig.colorbar(pc,ax=ax,shrink=0.5)
        
        
        ax = axs.flat[1]
        ax.set_aspect('equal')
        pc = ax.pcolormesh(*km(*grid_bike.axes),smb_bike,cmap='bwr_r',vmin=-1,vmax=1)
        fig.colorbar(pc,ax=ax,shrink=0.5)
        cf(ax, lon_bike, 360+np.arange(-60,-15,5))
        cf(ax, lat_bike, np.arange(60,90,5))
        cf(ax, frac_bike, [0.05], label=False)

        ax = axs.flat[2]
        ax.set_aspect('equal')
        pc = ax.pcolormesh(*km(*grid_bike.axes),sfct_bike,cmap='plasma',vmin=-50,vmax=0)
        fig.colorbar(pc,ax=ax,shrink=0.5)
        cf(ax, lon_bike, 360+np.arange(-60,-15,5))
        cf(ax, lat_bike, np.arange(60,90,5))
        cf(ax, frac_bike, [0.05], label=False)
        
        ax = axs.flat[3]
        fr =  np.sum(np.where(np.isfinite(frac_um),frac_um,0),axis=0)
        pc = ax.pcolormesh(*grid_um.axes,fr,
                           cmap='tab20c',vmin=0,vmax=1)
        ax.set_title(f'{np.min(fr)} <= frac col. tot. <= {np.max(fr)}')
        fig.colorbar(pc,ax=ax,shrink=0.5)
        ax.set_ylim(55,85)
        ax.set_xlim(260,355)