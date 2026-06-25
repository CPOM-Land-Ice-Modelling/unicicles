#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 15:50:48 2025

@author: ggslc
"""

import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator


def interp_to_local(zg, global_coords, local_coords, order=1, valid=None):
    """Interpolate 3D field from global to local grid coordinates.
    
    Interpolates a 3D field defined on a global grid to a local grid using
    either linear or nearest-neighbor interpolation. Handles missing data by
    filling null regions with nearest-neighbor values before interpolating.
    
    Parameters
    ----------
    zg : ndarray
        3D field on global grid, shape (nlev, N, M).
    global_coords : tuple of ndarray
        (xg, yg) coordinate arrays of global grid, each shape (N, M).
    local_coords : tuple of ndarray
        (xl, yl) coordinate arrays of local grid, each with same shape.
    order : int, optional
        Interpolation order: 1 for linear (default), 0 for nearest-neighbor.
    valid : ndarray, optional
        Boolean mask of valid data points on global grid, shape (nlev, N, M).
        If provided, null regions are filled with nearest-neighbor interpolation.
    
    Returns
    -------
    ndarray
        3D field interpolated to local grid, shape (nlev, *xl.shape).
    
    Raises
    ------
    ValueError
        If xl.shape != yl.shape or order not in (0, 1).
    """
    xg, yg = global_coords
    xl, yl = local_coords

    if xl.shape != yl.shape:
        raise ValueError('xl.shape != yl.shape')

    zl = np.empty((zg.shape[0],*xl.shape))

    if order not in (0,1):
        raise ValueError('order not in (0,1)')
        


    interpf = LinearNDInterpolator if order == 1 else NearestNDInterpolator

    for ec in range(0, zg.shape[0]):

        z = zg[ec,:,:]

        if isinstance(valid, np.ndarray):
            # fill the null regions with nearest-neigbours
            v = valid[ec,:,:]
            z = NearestNDInterpolator(
                np.array([xg[v].flat, yg[v].flat]).T, z[v].flat) (xg,yg)

        zl[ec,:,:] = interpf(np.array([xg.flat, yg.flat]).T, z.flat)(xl,yl)

    return zl

def interp_to_surface(field_l, topo_l, topo, lapse=0.0):
    """Interpolate field from pressure levels to surface topography.
    
    Vertically interpolates a field defined at pressure levels to the actual
    surface topography using linear interpolation between adjacent levels.
    Handles extrapolation above/below model levels using an optional lapse rate.
    
    Parameters
    ----------
    field_l : ndarray
        3D field at pressure levels, shape (nlev, N, M).
    topo_l : ndarray
        Topography (or altitude) at pressure levels, shape (nlev, N, M).
    topo : ndarray
        Target surface topography, shape (N, M).
    lapse : float, optional
        Lapse rate for extrapolation (K/m or appropriate units). Default is 0.0
        (no extrapolation adjustment).
    
    Returns
    -------
    ndarray
        Field interpolated to surface topography, shape (N, M).
    
    Notes
    -----
    - Levels are assumed ordered from lowest (0) to highest (-1) in index.
    - Linear interpolation is applied between adjacent levels where surface
      topography lies between their altitudes.
    - Extrapolation above the model top and below the model bottom uses
      the lapse rate to adjust values.
    """
    field = np.zeros(topo.shape)
    nec = topo_l.shape[0]

    for ec in range(1,nec):
        w = (topo - topo_l[ec,:,:])/(topo_l[ec-1,:,:] - topo_l[ec,:,:])
        wlo = np.where((topo >= topo_l[ec-1,:,:]) & (topo < topo_l[ec,:,:]),w,0)
        wup = np.where((topo >= topo_l[ec-1,:,:]) & (topo < topo_l[ec,:,:]),1-w,0)
        field += wup * field_l[ec,:,:] + wlo * field_l[ec-1,:,:]


    field = np.where(topo < topo_l[0,:,:],
                     field_l[0,:,:] + lapse*(topo_l[0,:,:] - topo),
                     field )

    field = np.where(topo > topo_l[-1,:,:],
                     field_l[-1,:,:] + lapse*(topo_l[-1,:,:] - topo),
                     field )


    return field
