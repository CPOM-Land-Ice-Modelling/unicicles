#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 10:14:16 2025

@author: ggslc
"""

import numpy as np
from pyproj import Proj
from scipy.interpolate import NearestNDInterpolator
from dataclass import dataclass
missing = np.nan
#widely used projections
PROJ_ANTARCTIC_3031 = '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
PROJ_ARCTIC_4326 = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'

def cell_id(i, j, ni, nj):
    return i + ni * j

class Box:
    """
    Box defined by lower and upper corner vectors lo, hi
    """

    def __init__(self, lo, hi):
        self._lo = lo
        self._hi = hi

    @property
    def lo(self):
        """
        lower corner vector
        """
        return self._lo

    @property
    def hi(self):
        """
        upper corner vector
        """
        return self._hi


class Uniform2DGrid:
    """
    2D grid with uniform spacing
    """

    def __init__(self, x, y):

        #check x, y ascending
        for xi in (x, y):
            if np.any(xi[:-1] >= xi[1:]):
                raise ValueError('x, y must be strictly ascending')

            dxi = xi[1] - xi[0]
            eps = 1.0e-10
            if np.any(xi[1:] > xi[:-1] + dxi + eps):
                raise ValueError('x, y must be uniformly spaced')


        #true for any grid defined by np.meshgrid(x, y)
        self._axes = (x, y)
        self._shape = (x.shape[0], y.shape[0])
        self._axes_index = tuple(np.arange(0,n) for n in self._shape)
        self._coords = np.meshgrid(x, y)
        self._array_shape = self._coords[0].shape

        # only for a uniform grid
        self._spacing = tuple(xi[1]-xi[0] for xi in self._axes)


    @property
    def array_shape(self):
        return self._array_shape

    @property
    def axes(self):
        return self._axes

    @property
    def axes_index(self):
        return self._axes_index

    @property
    def spacing(self):
        return self._spacing


    @property
    def coords(self):
        return self._coords

    def crop(self, crop_box):
        #axes = (x[np.where((x >= l) & (x <= h))] \
        #        for x,l,h in zip(self._axes, crop_box.lo, crop_box.hi))

        axes = []
        for i, x in enumerate(self._axes):
            t = x[np.where((x >= crop_box.lo[i]) & (x <= crop_box.hi[i]))]
            axes.append(t)
        x, y = axes[0], axes[1]
        return Uniform2DGrid(x, y)

def local_to_global_map_up(up_transform, global_grid, local_grid):

    xylg = up_transform(*local_grid.coords)
    ilg, jlg = [ np.intp((xlg-xg[0]) / dxg - 0.5) 
                for xlg, xg, dxg in zip(xylg, global_grid.axes, global_grid.spacing)]
    ng, mg = global_grid.array_shape
    return cell_id(ilg, jlg, mg, ng)

def local_to_global_map_down(down_transform, global_grid, local_grid):

    ig, jg = global_grid.axes_index
    ng, mg = global_grid.array_shape
    ii, jj = np.meshgrid(ig, jg)
    ij = cell_id(ii, jj, mg, ng)
    xyg = down_transform(*global_grid.coords)
    xyl = local_grid.coords
    return np.intp(NearestNDInterpolator(
        (np.array([xyg[0].flat, xyg[1].flat]).T), ij.flat)(*xyl))

def local_to_global_map(up_transform, down_transform, 
                        global_grid, local_grid, method="down"):
    
    if method == "down":
        return local_to_global_map_down(down_transform, 
                                        global_grid, local_grid)
    elif  method == "up":
        return local_to_global_map_up(up_transform, 
                                       global_grid, local_grid)
    else:
        raise ValueError(f'unknown local_to_global_map method {method}')

    return None


class GlobalLocal:
    
    def __init__(self, global_grid, local_grid, up_transform, down_transform):
        self._global_grid = global_grid
        self._local_grid = local_grid
        self._up_transform = up_transform
        self._down_transform = down_transform
        self._local_to_global_map = local_to_global_map(up_transform, down_transform, 
                        global_grid, local_grid)



def up_down_pair(proj):

    if isinstance(proj, str):
        if proj[0:5] == "+proj":
            #seems to be a proj string, so
            proj = Proj(proj)
        else:
            raise ValueError(f'unrecognised proj string {proj}')

    if isinstance(proj, Proj):

        def lon360(lon, lat):
            return lon%360, lat

        def up(x, y):
            return  lon360(*proj(x, y, inverse=True))

        def down(lon, lat):
            return proj(lon, lat, inverse=False)
    else:
        raise TypeError(f'unrecogised projection type {type(proj)}')

    return up, down


def grown_grid(downscale_transform, global_grid, local_grid):
    """
    Define a local (x,y) grid which covers at least the 
    local regions covered the supplied transformed global grid,
    and has cell centers and spacing in common with
    the supplied local (x,y) grid

    Parameters
    ----------
    downscale_transform : function(lon, lat) -> (x, y)
       tranformation mapping global (lon, lat) to local (x,y) co-ordinates
    global_grid : Uniform2DGrid
        global (lon, lat) grid spec - need not cover the globe 
    local_grid : TYPE
       local (x, y) grid spec

    Returns
    -------
    UniformGrid2D
        
    """

    xyg = downscale_transform(*global_grid.coords)
    xyl = local_grid.axes
    dxy = local_grid.spacing

    m = [max(0, int((- np.min(xg) + np.min(xl))/dxl) + 1)
         for xg, xl, dxl in zip(xyg, xyl, dxy)]

    n = [max(0, int((np.max(xg) - np.max(xl))/dxl) + 1)
         for xg, xl, dxl in zip(xyg, xyl, dxy)]

    p, q = [np.arange(x[0] - mx*dx, x[-1] + nx*dx, step=dx)
            for x, mx, nx, dx in zip(xyl, m, n, dxy)]

    return Uniform2DGrid(p, q)


def fraction_covered(up_transform, down_transform, global_grid, local_grid):
    """
    
    Compute fraction of each global grid cell covered by local grid cells

    Parameters
    ----------
    up_transform : function(x, y) -> (lon, lat)
        tranformation mapping local (x,y) to global (lon, lat) co-ordinates
    down_transform : function(lon, lat) -> (x, y)
        tranformation mapping global (lon, lat) to local (x,y) co-ordinates
    global_grid : Uniform2DGrid
        global (lon, lat) grid spec - need not cover the globe 
    local_grid : Uniform2DGrid
        local (x, y) grid spec

    Raises
    ------
    TypeError
        if  global_grid, local_grid are not Uniform2DGrid

    Returns
    -------
    numpy.ndarray(dtype = numpy.float64)
        fractional coverage

    """

    if not isinstance(global_grid, Uniform2DGrid):
        raise TypeError('global_grid not a Uniform2DGrid')

    if not isinstance(local_grid, Uniform2DGrid):
        raise TypeError('local_grid not a Uniform2DGrid')

    map_a = local_to_global_map(up_transform, down_transform, global_grid, local_grid)

    map_b = local_to_global_map(up_transform, down_transform, global_grid,
                                grown_grid(down_transform,
                                           global_grid, local_grid))

    ng, mg = global_grid.array_shape

    def cellid(i,j):
        return cell_id(i,j,mg,ng)

    #frac_coverage = np.zeros(global_grid.array_shape)
    ig, jg = global_grid.axes_index
    indx = cellid(*np.meshgrid(ig, jg)).flat[:]
    a = np.bincount(map_a.flat, minlength=1 + np.max(indx))
    b = np.bincount(map_b.flat, minlength=1 + np.max(indx))

    frac_coverage = np.where(b[indx] > 0,
                             a[indx]/np.float64(b[indx]),
                             0.0).reshape(global_grid.array_shape)

    return frac_coverage
