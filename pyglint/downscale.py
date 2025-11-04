#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 15:50:48 2025

@author: ggslc
"""

import numpy as np
from scipy.interpolate import LinearNDInterpolator

def interp_to_local(q_field_g, xg, yg, xl, yl):

    nec = q_field_g.shape[2] 
    nyl, nxl = xl.shape 
    if xl.shape != yl.shape:
        raise('xl.shape != yl.shape')
    q_field_l = np.empty((nyl,nxl,nec))
    for ec in range(0, nec):
        interp = LinearNDInterpolator(np.array([xg.flat, yg.flat]).T, 
                                      q_field_g[:,:,ec].flat)
        q_field_l[:,:,ec] = interp(xl,yl)
    return q_field_l

def interp_to_surface(q_field_l, topo_l, topo, lapse=0.0):
    
    field = np.zeros(topo.shape)
    nec = topo_l.shape[2]

    for ec in range(1,nec):
        w = (topo - topo_l[:,:,ec])/(topo_l[:,:,ec-1] - topo_l[:,:,ec])
        wlo = np.where((topo >= topo_l[:,:,ec-1]) & (topo < topo_l[:,:,ec]),w,0)
        wup = np.where((topo >= topo_l[:,:,ec-1]) & (topo < topo_l[:,:,ec]),1-w,0)
        field += wup * q_field_l[:,:,ec] + wlo * q_field_l[:,:,ec-1]
        
    
    field = np.where(topo < topo_l[:,:,0], 
                     q_field_l[:,:,0] + lapse*(topo_l[:,:,0] - topo), 
                     field )

    field = np.where(topo > topo_l[:,:,-1], 
                     q_field_l[:,:,-1] - lapse*(topo_l[:,:,-1] - topo), 
                     field )
    
        
    return field