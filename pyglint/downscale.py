#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 15:50:48 2025

@author: ggslc
"""

import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator


def interp_to_local(zg, global_coords, local_coords, order=1, valid=None):

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
