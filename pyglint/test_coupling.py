#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 09:48:10 2025

@author: ggslc
"""

import numpy as np
from coupling import atm_to_ism, ism_to_atm
from transformation import up_down_pair, Uniform2DGrid
from transformation import PROJ_ARCTIC_4326, PROJ_ANTARCTIC_3031
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from time import time


def read_ism_nc(nc_file):
    nc_bike_in = Dataset(nc_file,'r')
    topo = nc_bike_in['orog'][:,:]
    thk = nc_bike_in['lithk'][:,:]
    frac = nc_bike_in['sftgif'][:,:]
    calv = nc_bike_in['licalvf'][:,:]
    grid = Uniform2DGrid(nc_bike_in['x'][:].data, 
                         nc_bike_in['y'][:].data)
    mask = np.where(frac > 1e-4, True, False)
    
    return grid, topo, thk, calv, frac, mask

def read_atm_nc(nc_file):
    
    nc_um_in = Dataset(nc_file,'r')
    ilo, ihi = 148, 188
    jlo, jhi = 119, 140
    ilo, ihi = 0, -1
    jlo, jhi = 0, -1
    grid = Uniform2DGrid(nc_um_in['longitude'][ilo:ihi].data,
                            nc_um_in['latitude'][jlo:jhi].data)

    def prep_um(arr):
        return arr[:,jlo:jhi,ilo:ihi]
    
    area = prep_um(nc_um_in['tile_surface_area'])
    spy = 3600*360*24
    smb = prep_um(nc_um_in['ice_smb']) / 918 * spy # values suggest kg/s
    stemp = prep_um(nc_um_in['ice_stemp']) + 273.15
    snow = prep_um(nc_um_in['nonice_snowdepth'])
    shflx = prep_um(nc_um_in['snow_ice_hflux'])
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
        
    return grid, area, topo, topo_max, smb, stemp, snow, shflx
            

def cf(ax, x, y, z, zl, label=True):
    cs = ax.contour(x, y, z, zl, colors=['k'] ,linewidths=0.5, linestyles='-')
    if label:
        ax.clabel(cs, cs.levels, fontsize=8)

def km(x,y):
    return (x*1.0e-3,y*1e-3)


def plot_atm_coupler_inputs(grid, area, topo, topo_max, 
                            smb, stemp, snow, shflx):
    
    fig, axs = plt.subplots(4,1,figsize=(8,11), sharex=True)
    fig.suptitle('Coupler inputs <- atmosphere')
    
    ax = axs.flat[0]
    pc = ax.pcolormesh(*grid.axes, np.mean(stemp, axis=0),
                       vmin=230, vmax=275, cmap='plasma')
    fig.colorbar(pc,ax=ax,shrink=0.5, label=r'$T$ (K)')
    ax.set_title('Vertical mean ice surface temp.')
    
    ax = axs.flat[1]
    pc = ax.pcolormesh(*grid.axes, np.sum(smb*area, axis=0)/np.sum(area, axis=0), 
                       cmap='bwr_r', vmin=-1, vmax=1, )
    fig.colorbar(pc,ax=ax,shrink=0.5, label = r'SMB (m/a)')
    ax.set_title('Vertically integrated SMB')
    
    ax = axs.flat[2]
    pc = ax.pcolormesh(*grid.axes, np.sum(shflx*area, axis=0)/np.sum(area, axis=0), 
                       cmap='bwr_r', vmin=-.1, vmax=.1, )
    fig.colorbar(pc,ax=ax,shrink=0.5, label = r'$F_z$ (?)')
    ax.set_title('Vertically integrated heat flux')
    
    ax = axs.flat[3]
    pc = ax.pcolormesh(*grid.axes, np.sum(snow*area, axis=0)/np.sum(area, axis=0), 
                       cmap='bwr_r', vmin=-50, vmax=50, )
    fig.colorbar(pc,ax=ax,shrink=0.5, label = r'$h_s$ (m)')
    ax.set_title('Vertically integrated snow depth')
    
def plot_ism_heatmap(fig, ax, grid, lon, lat, var,  
                     title='', cmap='bwr',vmin=None, vmax=None,):
    
    ax.set_aspect('equal')
    ax.set_yticks([])
    ax.set_title(title)
    pc = ax.pcolormesh(*km(*grid.axes),var,cmap=cmap,
                       vmin=vmin,vmax=vmax)
    fig.colorbar(pc,ax=ax,shrink=0.5)
    #lon_lab = np.arange(lon.min().round(-1), lon.max().round(-1),30)
    #cf(ax,*km(*grid.axes), lon, lon_lab)
    #lat_lab = np.arange(lat.min().round(-1), lat.max().round(-1),10)
    #cf(ax,*km(*grid.axes),lat, lat_lab)
    #cf(ax, frac, [0.05], label=False)    
    
    
def plot_ism_coupler_inputs(grid, topo, thk, calv, frac, mask, up_tr):
    ... 
    
    def m(z):
        return np.ma.masked_array(z, frac < 1.0e-10)
    
    fig, axs = plt.subplots(2,2, figsize=(8,11))
    fig.suptitle('Coupler inputs <- ice sheet ')
    lon, lat = up_tr(*grid.coords)
    
    
    
    plot_ism_heatmap(fig, axs.flat[0], grid, lon, lat, m(topo), title='topo',
                 vmin=0,vmax=5000, cmap='plasma')
    plot_ism_heatmap(fig, axs.flat[1], grid, lon, lat, m(calv), 
                     title='calving flux',cmap='bwr',vmin=-1,vmax=1)
    plot_ism_heatmap(fig, axs.flat[2], grid, lon, lat, m(frac), 
                 title=r'ice frac',vmin=0,vmax=1, cmap='tab20c')
    plot_ism_heatmap(fig, axs.flat[3], grid, lon, lat, m(thk), 
                 title=r'ice thickness',vmin=0,vmax=5000, cmap='plasma')

    



def plot_ism_coupler_outputs(grid, smb, stemp, delta_snow, shflx, up_tr):
    
    lon, lat = up_tr(*grid.coords)
    
    fig, axs = plt.subplots(2,2, figsize=(8,11))
    fig.suptitle('Coupler outputs -> ice sheet ')
    plot_ism_heatmap(fig, axs.flat[0], grid, lon, lat, smb, title='SMB',
                     vmin=-1,vmax=1, cmap='bwr_r')
    plot_ism_heatmap(fig, axs.flat[1], grid, lon, lat, stemp, title='Temp',
                     vmin=220,vmax=275, cmap='bwr')
    plot_ism_heatmap(fig, axs.flat[2], grid, lon, lat, delta_snow, 
                     title=r'$\Delta$ snow',vmin=-50,vmax=50, cmap='bwr_r')
    plot_ism_heatmap(fig, axs.flat[3], grid, lon, lat, shflx, 
                     title=r'heat flux',vmin=-.1,vmax=.1, cmap='bwr_r')



    
def plot_atm_coupler_outputs(grid, area,  frac, topo, delta_snow, hflx, calv):    
    fig, axs = plt.subplots(3,1,figsize=(8,11), sharex=True)
    fig.suptitle('Coupler outputs -> atmosphere')
    
    ax = axs.flat[0]
    fr =  np.sum(np.where(np.isfinite(frac),frac,0),axis=0)
    pc = ax.pcolormesh(*grid.axes,fr,
                       cmap='tab20c',vmin=0,vmax=1)
    ax.set_title(f'{np.min(fr):2.2f} <= frac col. tot. <= {np.max(fr):2.2f}')
    fig.colorbar(pc,ax=ax,shrink=0.5)
    
    vi_area = np.sum(np.where(np.isfinite(area),area,0),axis=0) + 1.0e-10
    
    ax = axs.flat[1]
    cr =  np.sum(np.where(np.isfinite(calv),calv*area,0),axis=0)/vi_area
    pc = ax.pcolormesh(*grid.axes,cr,
                       cmap='Reds',vmin = 0, vmax=1e3)
    ax.set_title(f'{np.min(cr):2.2f} <= calving flux <= {np.max(cr):2.2f}')
    fig.colorbar(pc,ax=ax,shrink=0.5)
    
    ax = axs.flat[2]
    ds =  np.sum(np.where(np.isfinite(delta_snow),delta_snow*area,0),axis=0)/vi_area
    pc = ax.pcolormesh(*grid.axes,cr,
                       cmap='bwr',vmin = -50, vmax=50)
    ax.set_title(f'{np.min(ds):2.2f} <= delta snow <= {np.max(ds):2.2f}')
    fig.colorbar(pc,ax=ax,shrink=0.5)

#%%

plot = True

#%%
#transformation pair
up_tr_gris, down_tr_gris = up_down_pair(PROJ_ARCTIC_4326)

#example ice sheet model output (input to the coupler). bisicles CF
#Greenland
grid_gris, topo_gris, thk_gris, calv_gris, frac_gris, mask_gris = \
    read_ism_nc('bike_cf_gris4326.nc')

if plot:
    plot_ism_coupler_inputs(grid_gris, topo_gris, thk_gris, 
                            calv_gris, frac_gris, mask_gris, up_tr_gris)
    plt.savefig('bike_to_coupler_gris.png')

#Antarctoca
up_tr_ais, down_tr_ais = up_down_pair(PROJ_ANTARCTIC_3031)
grid_ais, topo_ais, thk_ais, calv_ais, frac_ais, mask_ais = \
    read_ism_nc('bike_cf_ais3031.nc')

if plot:
    plot_ism_coupler_inputs(grid_ais, topo_ais, thk_ais, 
                            calv_ais, frac_ais, mask_ais, up_tr_ais)
    plt.savefig('bike_to_coupler_ais.png')

#example atmosphere model output (input to the coupler).
grid_um, area_um, topo_um, topo_max, smb_um, stemp_um, snow_um, shflx_um,  = \
    read_atm_nc('atmos_cx209c_P1Y_20000101-20010101_icecouple.nc')
calv_um = np.ma.masked_array(np.zeros(smb_um.shape),smb_um.mask)

if plot:
    plot_atm_coupler_inputs(grid_um, area_um, topo_um, topo_max, 
                                smb_um, stemp_um, snow_um, shflx_um)

    plt.savefig('um_to_coupler.png')
#%% downscaling & upscaling

t = time()

delta_t = 1.0 # typical?

smb_gris, stemp_gris, delta_snow_gris, shflx_gris \
    = atm_to_ism(smb_um, stemp_um, snow_um, shflx_um,
                 topo_um, area_um, grid_um,
                 topo_gris, thk_gris, frac_gris, mask_gris, grid_gris,
                 up_tr_gris, down_tr_gris, delta_t)

print(time() - t)

if plot:
    plot_ism_coupler_outputs(grid_gris, smb_gris, stemp_gris, 
                             delta_snow_gris, shflx_gris, up_tr_gris)
    plt.savefig('coupler_to_bike_gris.png')

t = time()

smb_ais, stemp_ais, delta_snow_ais, shflx_ais \
    = atm_to_ism(smb_um, stemp_um, snow_um, shflx_um,
                 topo_um, area_um, grid_um,
                 topo_ais, thk_ais, frac_ais, mask_ais, grid_ais,
                 up_tr_ais, down_tr_ais, delta_t)
    
print(time() - t)

if plot:
    plot_ism_coupler_outputs(grid_ais, smb_ais, stemp_ais, 
                             delta_snow_ais, shflx_ais, up_tr_ais)
    plt.savefig('coupler_to_bike_ais.png')


um_shape = smb_um.shape
frac_um = np.full(um_shape, 0.0)

t = time()


frac_um, topo_um, delta_snow_um, hlfx_um, calv_um, \
    = ism_to_atm(topo_max, area_um, grid_um,  
                     frac_um, snow_um, shflx_um, calv_um, topo_um,
                     topo_gris, frac_gris, 
                     shflx_gris, delta_snow_gris, calv_gris,
                     mask_gris, grid_gris,
                     up_tr_gris, down_tr_gris)

print(time() - t)


t = time()

frac_um, topo_um, delta_snow_um, hlfx_um, calv_um, \
    = ism_to_atm(topo_max, area_um, grid_um,  
                     frac_um, snow_um, shflx_um, calv_um, topo_um,
                     topo_ais, frac_ais, 
                     shflx_ais, delta_snow_ais, calv_ais,
                     mask_ais, grid_ais,
                     up_tr_ais, down_tr_ais)

print(time() - t)


if plot:
    plot_atm_coupler_outputs(grid_um, area_um, frac_um, topo_um, delta_snow_um, 
                             hlfx_um, calv_um)
    plt.savefig('coupler_to_um.png')