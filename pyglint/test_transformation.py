#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 13:05:34 2025

@author: ggslc
"""

from transformation import *
import matplotlib.pyplot as plt
import unittest

def map_and_fractions():
    
   
    up, down = up_down_pair(PROJ_ANTARCTIC_3031)

    km = 1.0e3
    nx = 256
    dx = 4000*km/nx
    x = np.linspace(-2000*km + dx/2, 2000*km - dx/2, nx)
    y = np.linspace(-2000*km + dx/2, 2000*km - 3*dx/2, nx-1)
    local_grid = Uniform2DGrid(x, y)


    nlon, nlat = 36, 15
    dlon, dlat = 360.0/nlon, 30/nlat
    lon = np.linspace(dlon/2.0, 360.0-dlon/2.0, nlon)
    lat = np.linspace(-90+dlat/2, -60-dlat/2, nlat)
    global_grid = Uniform2DGrid(lon, lat)

    lgm =  local_to_global_map(down, global_grid, local_grid)

    fr = fraction_covered(down, global_grid, local_grid)


    fig, axs = plt.subplots(1, 2, figsize=(8,4))

    ax= axs[0]
    ax.set_aspect('equal')
    ax.pcolormesh(x, y, lgm%20, cmap='tab20c')

    lon_ism, lat_ism = up(*local_grid.coords)

    def cf(ax, xx, yy, z, zl, label=True):
        cs = ax.contour(xx, yy, z, zl,
                       colors=['k'] ,linewidths=0.5, linestyles='-')
        if label:
            ax.clabel(cs, cs.levels, fontsize=8)

    cf(ax,x,y,lon_ism, lon + dlon/2, label=False  )
    cf(ax,x,y,lat_ism, lat + dlat/2,  label=False )



    ax = axs[1]
    pc = ax.pcolormesh(lon, lat, fr, vmin=0, vmax=2, cmap='bwr')
    fig.colorbar(pc, ax = ax)
    ax.set_yticks(lat + dlat/2)
    ax.set_xticks(lon + dlon/2)
    ax.grid(color='k',lw=0.5)
    
    return lgm, fr

class TestTranformation(unittest.TestCase):    

    def test_fractions(self):
         self.assertAlmostEqual(np.max(arr_frac),1.0)
         self.assertAlmostEqual(np.min(arr_frac),0.0)
         self.assertAlmostEqual(np.mean(arr_frac), 0.6783907677463454)
         
    def test_map(self):
         self.assertEqual(np.min(arr_map), 0)
         self.assertEqual(np.max(arr_map), 463)
         self.assertEqual(np.unique(arr_map).shape, (396,))
         
if __name__ == "__main__":

   arr_map, arr_frac =  map_and_fractions()
   unittest.main()
    