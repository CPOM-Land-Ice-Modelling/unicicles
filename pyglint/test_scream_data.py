#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 13:01:02 2025

@author: ggslc
"""

import numpy as np
import unittest
import matplotlib.pyplot as plt
from scipy.interpolate import NearestNDInterpolator
from transformation import cell_id

LAPSE_RATE = 0.0085

def scream_surface_test_data(*args, local_grid_shape = (511,512), 
                             global_grid_shape = (32,31),
                             global_min_elev = 0.0,
                             global_max_elev = 1500,
                             global_grid_angle = np.pi/6.0,
                             n_elev = 5, 
                             lapse=LAPSE_RATE, 
                             z_freeze = 1000.0, 
                             z_freeze_var = 250.0,
                             usrf_min = 0.0, usrf_max = 2000.0, 
                             **kwargs ):

    # test data, elevation classes
    topo_max_ec = np.linspace(global_min_elev, global_max_elev, n_elev+1)
    
    # test data / local grid. centers.
    dx = 1.0
    m, n = local_grid_shape
    y = np.linspace(-dx*m/2+dx/2, dx*m/2-dx/2, m)
    x = np.linspace(-dx*n/2+dx/2, dx*n/2-dx/2, n)
    xx, yy = np.meshgrid(x, y)
    
    # test surface & mask on local grid
    usrf = 0.7*np.cos(14*xx/m) + 0.3*np.sin(12*yy/n)
    mask = np.where(xx**2+yy**2 < (.45*min(m,n))**2,True,False)
    usrf = np.where(mask, usrf_min + usrf_max * np.abs(usrf), 0.0)
    xxx, yyy, zzz = np.meshgrid(x, y, topo_max_ec[:-1])
    
    
    #test data / global 2D grid
    M, N = global_grid_shape
    DX = m/M
    Y = np.linspace(-DX*M*3/4, DX*M*3/4, M)
    X = np.linspace(-DX*N*3/4, DX*N*3/4, N)
    XX, YY = np.meshgrid(X, Y)
   
   
    
    
    #local co-ordinates of global grid points
    alpha = global_grid_angle
    xXY = XX * np.cos(alpha) - YY * np.sin(alpha)
    yXY = XX * np.sin(alpha) + YY * np.cos(alpha)
    
    
    #global co-ordinates of local grid points
    Xxy = xx * np.cos(alpha) + yy * np.sin(alpha)
    Yxy = -xx * np.sin(alpha) + yy * np.cos(alpha)
    
    #local to global index map - a cell id for every local grid cell
    J = np.arange(0,M)
    I = np.arange(0,N)
    II, JJ = np.meshgrid(I, J)
    CC = cell_id(II, JJ, N, M)
    local_to_global_map = NearestNDInterpolator((np.array([xXY.flat, yXY.flat]).T), CC.flat)(xx,yy)

    #test sftc_g
    sftc_g = np.zeros([n_elev, M, N])
    topo_g = np.zeros([n_elev, M, N])
    
    L = np.max(x) - np.min(x)
    #Z = np.cos(8*xXY/L) + np.sin(9*yXY/L)
    #z = np.cos(8*xx/L) + np.sin(9*yx/L)
    #s_freeze = 1000.0# + 250.0*Z
    #lapse = 9.0e-3
    zfun = lambda x,y :  z_freeze_var*(np.cos(8*x/L) + np.sin(9*y/L))
    Tfun = lambda s, x, y: -lapse*(s-z_freeze + zfun(x,y)) 
    
    for ec in range(0, n_elev):
        topo_g[ec,:,:] = 0.5*(topo_max_ec[ec]+topo_max_ec[ec+1])
        sftc_g[ec,:,:] = Tfun(topo_g[ec,:,:],xXY,yXY)
    sftc_s =  Tfun(usrf,xx,yy) 
        
    return  X, Y, XX, YY, xXY, yXY, x, y, xx, yy, \
        local_to_global_map, sftc_g, topo_g, usrf, sftc_s, \
        mask,topo_max_ec 