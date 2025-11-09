#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 10:14:16 2025

@author: ggslc
"""

import numpy as np
missing = np.nan

from pyproj import Proj

#widely used projections
PROJ_ANTARCTIC_3031 = '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
PROJ_ARCTIC_4326 = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'

def cell_id(i, j, ni, nj):
    return j + nj * i

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


def UpDown(proj):
    
    if type(proj) is str:
        if proj[0:5] == "+proj":
            #seems to be a proj string, so 
            proj = Proj(proj)
    
    if type(proj) is Proj:
        up = Upscale(lambda x,y : proj(x,y,inverse=True))
        down = Downscale(lambda x,y : proj(x,y,inverse=False))
    else:
        raise ('unrecogised projection')
        
    return up, down


if __name__ == "__main__":
 
    import matplotlib.pyplot as plt   
    up, down = UpDown(PROJ_ANTARCTIC_3031)
    
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