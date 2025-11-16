#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 10:14:16 2025

@author: ggslc
"""

import numpy as np
from pyproj import Proj
from scipy.interpolate import NearestNDInterpolator

missing = np.nan
#widely used projections
PROJ_ANTARCTIC_3031 = '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
PROJ_ARCTIC_4326 = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'

def cell_id(i, j, ni, nj):
    return j + nj * i

class Box:

    def __init__(self, lo, hi):
        self._lo = lo
        self._hi = hi

    @property
    def lo(self):
        return self._lo

    @property
    def hi(self):
        return self._hi


class Uniform2DGrid:


    def __init__(self, x, y):

        self._axes = (x, y)
        self._shape = tuple(len(xi) for xi in self._axes)
        self._axes_index = tuple(np.arange(0,N) for N in self._shape)
        self._coords = np.meshgrid(x, y)

    @property
    def shape(self):
        return self._shape

    @property
    def axes(self):
        return self._axes

    @property
    def axes_index(self):
        return self._axes_index


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

        return Uniform2DGrid(*axes)


class Downscale:

    def __init__(self, transform):
        self._transform = transform

    def local_xy(self, global_XY):
        return self._transform(*global_XY)

class Upscale:

    def __init__(self, transform):
        self._transform = transform

    def global_XY(self, local_xy):
        return self._transform(*local_xy)


def local_to_global_map(downscale_transform, global_grid, local_grid):

    I, J = global_grid.axes_index
    N, M = global_grid.shape
    II, JJ = np.meshgrid(I, J)
    CC = cell_id(II, JJ, N, M)
    XY = downscale_transform.local_xy(global_grid.coords)
    xy = local_grid.coords
    return NearestNDInterpolator(
        (np.array([XY[0].flat, XY[1].flat]).T), CC.flat)(*xy)

def lon_fiddle(lonlat):
    return lonlat[0]%360,lonlat[1]

def up_down_pair(proj):

    if isinstance(proj, str):
        if proj[0:5] == "+proj":
            #seems to be a proj string, so
            proj = Proj(proj)
        else:
            raise ValueError(f'unrecognised proj string {proj}')

    if isinstance(proj, Proj):
        post_proj = lambda X,Y : (X%360, Y)
        up = Upscale(lambda x,y : post_proj(*proj(x,y,inverse=True)))
        down = Downscale(lambda x,y : proj(x,y,inverse=False))
    else:
        raise TypeError(f'unrecogised projection type {type(proj)}')

    return up, down


if __name__ == "__main__":

    import matplotlib.pyplot as plt
    up, down = up_down_pair(PROJ_ANTARCTIC_3031)

    km = 1.0e3
    x = np.linspace(-3000*km, 3000*km, 81)
    y = np.linspace(-3000*km, 3000*km, 81)

    xx, yy =  np.meshgrid(x, y)
    llon,llat = up.global_XY((xx, yy))

    xxx, yyy = down.local_xy((llon, llat))

    z = xx**2 + yy**2

    plt.pcolormesh(x,y,z)
    plt.contour(x,y,xxx)
    plt.contour(x,y,yyy)
    