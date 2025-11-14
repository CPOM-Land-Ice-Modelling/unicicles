#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 15:50:48 2025

@author: ggslc
"""

import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.interpolate import RegularGridInterpolator

def interp_to_local(field_g, global_coords, local_coords, order=1, valid=None):

    xg, yg = global_coords
    xl, yl = local_coords

    nec = field_g.shape[0] 
    nyl, nxl = xl.shape 
    if xl.shape != yl.shape:
        raise('xl.shape != yl.shape')
    field_l = np.empty((nec,nyl,nxl))
    
    if not (order == 1 or order == 0):
        raise('order != 0 or 1')
    
    interpf = LinearNDInterpolator if order == 1 else NearestNDInterpolator
    
    for ec in range(0, nec):

        z =  field_g[ec,:,:]      

        if type(valid) is np.ndarray:
            # fill the null regions with nearest-neigbours
            v = valid[ec,:,:]
            zf = NearestNDInterpolator(
                np.array([xg[v].flat, yg[v].flat]).T, z[v].flat)
            z = zf(xg,yg)
                
        field_l[ec,:,:] = interpf(np.array([xg.flat, yg.flat]).T, z.flat)(xl,yl)
        
        
        
    return field_l

def interp_to_surface(field_l, topo_l, topo, lapse=0.0):
    
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