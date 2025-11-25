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

def glint_conservation_adjust(f2d_ism_nc, mask_ism,  
                              f3d_atm, area_atm,  valid_atm,
                              grid_pair):



    # conservation enforcement - see glint_downscaling_gcm
    # roughly - adjust smb so that the total smb_ism within an
    # atmosphere grid box equals the smb_atm integrated vertically.
    # this resembles the glint code. blocky.
    # there is a NASTY! comment in the glint code
    f2d_ism = np.zeros(f2d_ism_nc.shape) + f2d_ism_nc

    ism_to_atm_map = grid_pair.local_to_global_map
    grid_ism = grid_pair.local_grid
    grid_atm = grid_pair.global_grid


    dx, dy = grid_ism.spacing # \todo - move to Uniform2DGrid

    smb_adjust = np.zeros(grid_atm.array_shape)
    N, M  = grid_atm.array_shape
    for I in grid_atm.axes_index[0]:
        for J in grid_atm.axes_index[1]:
            if np.any(valid_atm[:, J, I]):

                cc = cell_id(I, J, M, N)

                # vertical sum
                atm_col_total = np.sum(f3d_atm[:, J, I]*area_atm[:, J, I])

                # horizontal sum
                ism_total = dx**2 * np.sum(np.where(ism_to_atm_map == cc, f2d_ism, 0.0))
                ism_area =  dx**2 * np.sum(np.where(ism_to_atm_map == cc, 1.0, 0.0))

                if ism_area > 0:
                    smb_adjust[J, I] = (atm_col_total - ism_total)/ism_area
                    f2d_ism = np.where(ism_to_atm_map == cc,
                                       f2d_ism + smb_adjust[J, I],
                                       f2d_ism)

    return np.ma.masked_array(f2d_ism, ~mask_ism)


def crop_global(arrs_global, grid_atm, grid_ism, up_transform):

    # subset
    lon_ism, lat_ism = up_transform(*grid_ism.coords)
    dlon, dlat = ( np.max(np.abs(ll[1:] - ll[:-1])) for ll in  grid_atm.axes)
    bb = Box( [np.min(lon_ism) - dlon, np.min(lat_ism) - dlat],
             [np.max(lon_ism) + dlon, np.max(lat_ism) + dlat])


    ilo = np.argmin( np.abs(grid_atm.axes[1] - bb.lo[1])) - 1
    ihi = np.argmin( np.abs(grid_atm.axes[1] - bb.hi[1])) + 1
    grid_atm_sub = Uniform2DGrid( grid_atm.axes[0], grid_atm.axes[1][ilo:ihi])

    arrs = [arr[:,ilo:ihi,:] for arr in arrs_global]
    arrs.append(grid_atm_sub)
    arrs.append(ilo)
    arrs.append(ihi)
    return arrs





def atm_to_ism(smb_atm, sfct_atm, topo_atm, area_atm,  grid_atm,
               topo_ism, frac_ism, mask_ism, grid_ism,
               up_transform,  down_transform):


    smb_atm, sfct_atm, topo_atm, area_atm, grid_atm, ilo, ihi = \
        crop_global((smb_atm, sfct_atm, topo_atm, area_atm),
                    grid_atm, grid_ism, up_transform)

    grid_pair = GlobalLocalGridPair(grid_atm, grid_ism, up_transform, down_transform)  

    atm_coords = down_transform(*grid_atm.coords)

    valid_atm = area_atm.data > 0.0
    smb_xyz = interp_to_local(smb_atm.data, atm_coords, grid_ism.coords,
                              order=1, valid=valid_atm & ~smb_atm.mask)

    sfct_xyz = interp_to_local(sfct_atm.data, atm_coords, grid_ism.coords,
                               order=1, valid=~sfct_atm.mask)

    # this is probably uniform in x,y , but ...
    topo_xyz = interp_to_local(topo_atm, atm_coords, grid_ism.coords, order=1)


    sfct_ism, smb_ism = \
        (np.ma.masked_array(interp_to_surface(f_xyz, topo_xyz, topo_ism), ~mask_ism)
         for f_xyz in (sfct_xyz, smb_xyz))

    if not isinstance(area_atm, type(None)):
        smb_ism = glint_conservation_adjust(
            smb_ism, mask_ism, 
            smb_atm, area_atm, valid_atm,
            grid_pair)

    return smb_ism, sfct_ism


def splice_global(global_arr, region_arr, frac_cover,  ilo, ihi):

    global_arr[:,ilo:ihi,:] = region_arr*frac_cover + (1.0-frac_cover)*global_arr[:,ilo:ihi,:]

    return global_arr

def ism_to_atm(topo_max_atm, area_atm, grid_atm,
               ice_frac_atm, topo_atm,
               topo_ism, frac_ism, mask_ism, grid_ism,
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

    # \todo check splice_fields_gcm

    dbg = 0

    area_atm_sub, grid_atm_sub, ilo, ihi = \
        crop_global([area_atm], grid_atm, grid_ism, up_transform)
        
    grid_pair = GlobalLocalGridPair(grid_atm_sub, grid_ism, up_transform, down_transform)    
        
    region_data = mean_to_global_mec([frac_ism.data, topo_ism.data],
                                                [True, False],
                                                topo_ism.data, mask_ism,
                                                topo_max_atm, 
                                                grid_pair.local_to_global_map,
                                                area_atm_sub.shape)

    frac_cover = fraction_covered(grid_pair)
    global_data = [ice_frac_atm, topo_atm]
    global_data = [splice_global(g, r , frac_cover, ilo, ihi)
                   for g,r in zip(global_data, region_data)]

    return global_data





