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

def scream_surface_test_data(local_grid_shape = (512,512), 
                             global_grid_shape = (32,32),
                             global_min_elev = 0.0,
                             global_grid_angle = np.pi/6.0,
                             n_elev = 6):

    # test data, elevation classes
    topo_max_ec = np.linspace(0, 2400, n_elev+1)
    
    # test data / local grid. centers.
    dx = 1.0
    m, n = 512, 512
    y = np.linspace(-dx*m/2+dx/2, dx*m/2-dx/2, m)
    x = np.linspace(-dx*n/2+dx/2, dx*n/2-dx/2, n)
    xx, yy = np.meshgrid(x, y)
    
    # test surface & mask on local grid
    usrf = 0.7*np.cos(14*xx/m) + 0.3*np.sin(12*yy/n)
    mask = np.where(xx**2+yy**2 < (m/3)**2,1,0)
    usrf = mask * 2000.0 *np.abs(usrf)
    xxx, yyy, zzz = np.meshgrid(x, y, topo_max_ec[:-1])
    
    
    #test data / global 2D grid
    M, N = 32, 28
    DX = m/M
    Y = np.linspace(-DX*M*3/4, DX*M*3/4, M)
    X = np.linspace(-DX*N*3/4, DX*N*3/4, N)
    XX, YY = np.meshgrid(X, Y)
    L = np.max(X) - np.min(X)
    Z = np.cos(8*XX/L) + np.sin(9*YY/L)
    
    #XXX,YYY,ZZZ  = np.meshgrid(X, Y, topo_max_ec[:-1])
    
    #local co-ordinates of global grid points
    alpha = global_grid_angle
    xXY = XX * np.cos(alpha) - YY * np.sin(alpha)
    yXY = XX * np.sin(alpha) + YY * np.cos(alpha)
    
    
    #global co-ordinates of local grid points
    Xxy = xx * np.cos(alpha) + yy * np.sin(alpha)
    Yxy = -xx * np.sin(alpha) + yy * np.cos(alpha)
    
    DX = 1
    J = np.arange(0,M)
    I = np.arange(0,N)
    II, JJ = np.meshgrid(I, J)
    CC = cell_id(II, JJ, M, N)
    local_to_global_map = NearestNDInterpolator((np.array([xXY.flat, yXY.flat]).T), CC.flat)(xx,yy)

    
    #test q_sftc_g
    q_sftc_g = np.zeros([M, N, n_elev])
    q_topo_g = np.zeros([M, N, n_elev])
    freeze = 1000.0 + 250.0*Z
    lapse = 9.0e-3
    for ec in range(0, n_elev):
        q_topo_g[:,:,ec] = (topo_max_ec[ec]) + global_min_elev,
        q_sftc_g[:,:,ec] = -lapse*( q_topo_g[:,:,ec] - freeze)
        
        
    return  X, Y, XX, YY, xXY, yXY, x, y, xx, yy, \
        local_to_global_map, q_sftc_g, q_topo_g, usrf, mask,topo_max_ec 