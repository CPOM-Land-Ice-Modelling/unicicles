#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 16:10:16 2025

@author: ggslc
"""
import numpy as np
from downscale import interp_to_local, interp_to_surface
from upscale import mean_to_global_mec
from transformation import local_to_global_map, cell_id, fraction_covered
from transformation import  Uniform2DGrid, Box, GlobalLocalGridPair


def delta_h_snow(h_ice, h_snow, thresh=50.0):
    
    return np.where( h_ice < thresh,  
                      np.where ( h_snow > thresh, - h_snow, h_ice), 0.0)
    
def glint_conservation_adjust(f2d_ism_nc, mask_ism,  
                              f3d_atm, area_atm,  valid_atm,
                              grid_pair):

    
    # conservation enforcement - see glint_downscaling_gcm
    # roughly - adjust smb so that the total smb_ism within an
    # atmosphere grid box equals the smb_atm integrated vertically.
    # this resembles the glint code. blocky.
    # there is a NASTY! comment in the glint code
    # optmized with the help of co-pilot  

    f2d_ism = np.array(f2d_ism_nc, copy=True)

    ism_to_atm_map = grid_pair.local_to_global_map
    grid_ism = grid_pair.local_grid
    grid_atm = grid_pair.global_grid

    dx, dy = grid_ism.spacing 
    area_factor = dx * dy

    # Number of atmosphere cells expected (cell_id produces indices in [0, N*M-1])
    N, M = grid_atm.array_shape
    atm_cell_count = N * M

    # Compute vertical-integrated atm column totals safely: zero out invalid vertical entries
    # f3d_atm and area_atm shapes: (lev, N, M) -> sum over axis=0 results in (N, M)
    atm_col_total_2d = np.sum(np.where(valid_atm, f3d_atm * area_atm, 0.0), axis=0)

    # Flatten to length N*M with C-order so that index cc = J*M + I matches cell_id(I,J,M,N)
    atm_col_total_flat = atm_col_total_2d.ravel(order='C')

    # Prepare flattened local mappings and local smb
    mapped_ids = ism_to_atm_map.ravel(order='C').astype(np.int64)
    f2d_flat = f2d_ism.ravel(order='C')
    m2d_flat = np.where(mask_ism.ravel(order='C'),1.0,0.0)

    # Sum f2d_ism per atm-cell and count number of local cells per atm-cell
    ism_total_flat = area_factor * np.bincount(mapped_ids, weights=f2d_flat*m2d_flat, minlength=atm_cell_count)
    ism_area_flat = area_factor * np.bincount(mapped_ids, minlength=atm_cell_count)
    ism_ice_area_flat = area_factor * np.bincount(mapped_ids, weights=m2d_flat, minlength=atm_cell_count)

    # Compute adjustments only for atm cells that have any valid vertical data and have non-zero ism area
    atm_valid_flat = np.any(valid_atm, axis=0).ravel(order='C')

    delta_flat = np.zeros(atm_cell_count, dtype=float)
    nonz = atm_valid_flat & (ism_ice_area_flat > 0.0)
    if np.any(nonz):
        delta_flat[nonz] = (atm_col_total_flat[nonz] *  \
                                ism_ice_area_flat[nonz] / ism_area_flat[nonz]    \
                                    - ism_total_flat[nonz]) / ism_ice_area_flat[nonz]

    # Apply adjustments to local grid: each local cell gets the delta of its mapped atm cell
    # For unmapped local cells (flat_map < 0) we do not apply any adjustment (delta_contrib = 0)
    delta_per_local = np.zeros_like(f2d_flat, dtype=float)
    # Only assign for mapped local cells to avoid indexing errors
    delta_per_local = delta_flat[mapped_ids]

    # Update flattened f2d and reshape back
    f2d_flat += delta_per_local
    f2d_ism = f2d_flat.reshape(f2d_ism.shape, order='C')

    return np.ma.masked_array(f2d_ism, ~mask_ism)



def crop_global(arrs_global, grid_atm, grid_ism, up_transform):

    # subset
    lon_ism, lat_ism = up_transform(*grid_ism.coords)
    dlon, dlat = ( np.max(np.abs(ll[1:] - ll[:-1])) for ll in  grid_atm.axes)
    bb = Box( [np.min(lon_ism) - dlon, np.min(lat_ism) - dlat],
             [np.max(lon_ism) + dlon, np.max(lat_ism) + dlat])


    ilo = max(0,np.argmin( np.abs(grid_atm.axes[1] - bb.lo[1])) - 1)
    ihi = min(np.argmin( np.abs(grid_atm.axes[1] - bb.hi[1])) + 1,grid_atm.axes[1].shape[0])
    grid_atm_sub = Uniform2DGrid( grid_atm.axes[0], grid_atm.axes[1][ilo:ihi])

    arrs = [arr[:,ilo:ihi,:] for arr in arrs_global]
    arrs.append(grid_atm_sub)
    arrs.append(ilo)
    arrs.append(ihi)
    return arrs





def atm_to_ism(smb_atm, stemp_atm, snow_atm, shflx_atm,  
               topo_atm, area_atm,  grid_atm,
               topo_ism, lithk_ism, frac_ism, mask_ism, grid_ism,
               up_transform,  down_transform, time_step):

    
    
    smb_atm, stemp_atm, snow_atm, shflx_atm, topo_atm, area_atm, grid_atm, \
    ilo, ihi = \
        crop_global((smb_atm, stemp_atm, snow_atm, shflx_atm, topo_atm, area_atm),
                    grid_atm, grid_ism, up_transform)

    grid_pair = GlobalLocalGridPair(grid_atm, grid_ism, up_transform, down_transform)  

    atm_coords = down_transform(*grid_atm.coords)

    valid_atm = area_atm.data > 0.0
    smb_xyz, snow_xyz, shflx_xyz  = [interp_to_local(f_atm.data, atm_coords, grid_ism.coords,
                              order=1, valid=valid_atm & ~f_atm.mask) 
                         for f_atm in [smb_atm, snow_atm, shflx_atm]]

    stemp_xyz = interp_to_local(stemp_atm.data, atm_coords, grid_ism.coords,
                               order=1, valid=~stemp_atm.mask)

    # this is probably uniform in x,y , but ...
    topo_xyz = interp_to_local(topo_atm, atm_coords, grid_ism.coords, order=1)


    stemp_ism, smb_ism, snow_ism, shflx_ism = \
        (np.ma.masked_array(interp_to_surface(f_xyz, topo_xyz, topo_ism), ~mask_ism)
         for f_xyz in (stemp_xyz, smb_xyz, snow_xyz, shflx_xyz))

    if True and not isinstance(area_atm, type(None)):
        smb_ism, snow_ism, shflx_ism = [glint_conservation_adjust(
                            f_ism, mask_ism, 
                            f_atm, area_atm, valid_atm,
                            grid_pair) 
                for f_ism, f_atm in zip([smb_ism, snow_ism, shflx_ism],
                                 [smb_atm, snow_atm, shflx_atm])]
        
        
    #snow to ice conversion - store *change in* snow depth in snow_ism    
    delta_snow_ism = np.ma.masked_array(delta_h_snow(lithk_ism, snow_ism), ~mask_ism)
    smb_ism -= delta_snow_ism / time_step 

    return smb_ism, stemp_ism, delta_snow_ism, shflx_ism


def splice_global(global_arr, region_arr, frac_cover,  ilo, ihi):

    global_arr[:,ilo:ihi,:] = region_arr*frac_cover + (1.0-frac_cover)*global_arr[:,ilo:ihi,:]

    return global_arr

def ism_to_atm(topo_max_atm, area_atm, grid_atm,
               ice_frac_atm, delta_snow_atm, hflux_atm, calv_atm, topo_atm, 
               topo_ism, frac_ism, 
               hflx_ism, d_snow_ism, calv_ism,
               mask_ism, grid_ism,
               up_transform, down_transform):

    # GLINT 3D ouputs : gtopo - ice sheet surface,
    #                   gfrac - ice covered fraction
    #                   grofi - ice run off (calving?)
    #                   grofl - liquid run off
    #                   ghflx - heat flux
    #                   gsdep - snow depth (in total / out anomaly)
    #                   gcalv = calving flux
    # GLINT 0D outputs : ice volume

    # wrappers/ukesm-ice_NETCDF/gl_mod.f90 calls glint to obtain
    # gsdep (in also), gfrac, ghflx, gcalv, glfrac, ice_volume
    # i.e ignore grofi , grofl

    dbg = 0

    area_atm_sub, grid_atm_sub, ilo, ihi = \
        crop_global([area_atm], grid_atm, grid_ism, up_transform)
        
    grid_pair = GlobalLocalGridPair(grid_atm_sub, grid_ism, up_transform, down_transform)    
        

    region_data = mean_to_global_mec([frac_ism.data, 
                                      topo_ism.data, 
                                      d_snow_ism.data, 
                                      hflx_ism.data, 
                                      calv_ism.data], 
                                     [True, False, False, False, False], 
                                     topo_ism.data, mask_ism,
                                     topo_max_atm, 
                                     grid_pair.local_to_global_map,
                                     area_atm_sub.shape)
   
    
    frac_cover = fraction_covered(grid_pair)
    global_data = [ice_frac_atm, topo_atm, delta_snow_atm, hflux_atm, calv_atm]
    global_data = [splice_global(g, r , frac_cover, ilo, ihi)
                   for g,r in zip(global_data, region_data)]

    return global_data





