#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 30 16:13:44 2025

@author: ggslc
"""

import numpy as np
from functools import partial
from transformation import cell_id, missing


def subset_ec(topo, topo_max, ec):
    return (topo < topo_max[ec+1]) & (topo >= topo_max[ec])


def local_to_global_cell_mec(field, topo, mask, topo_max, local_to_global_map,
                            lcolfrac, cell_indx):

    indx = mask & (local_to_global_map  == cell_indx)

    nec = len(topo_max) - 1  # 0,z1,z2,,,,zn
    field_sum = np.zeros(nec)
    field_count = np.zeros(nec)
    ts = topo[indx]
    tf = field[indx]
    field_count_col = len(tf) 
    for ec in range(0, nec):
        kndx = subset_ec(ts, topo_max, ec)
        t = tf[kndx]
        field_sum[ec] = np.sum(t)
        field_count[ec] = len(t) 

    return field_sum, field_count, field_count_col


def mean_to_global_mec(field, topo, mask, topo_max, local_to_global_map,
                       global_shape, lcolfrac=False, missing = missing):
    """
    
    Compute the means of 2D field on the global 3D grid.
    Cell values will be non-zero where the 2D surface elevation (topo) 
    intersects the cell. Cell vertical boundaries are defined
    by topo_max
    
    Parameters
    ----------
    field : numpy.ndarray / float 2D array
        field on the local grid 
    topo : numpy.ndarray / float 2D array
        surface elevation on the local grid
    topo_max : numpy.ndarray / float 1D array
        elevation class limits
    local_to_global_map : 2D array with
        global grid cell ID for each local grid cell
    global_shape : (int, int, int)
        shape of the global grid
    lcolfrac : bool, optional
        If true, mean computed over global column, if false mean 
        computef over the gobal cell. The default is False.

    Returns
    -------
    global_field : numpy.ndarray / float 3D array
        means of input field for each global cell

    """

    nJ, nI, __ = global_shape
    global_field = np.full(global_shape, missing) 

    for I in range(0, nI):
        for J in range(0, nJ):
            fsum, count, count_col  = local_to_global_cell_mec( \
                    field, topo, mask, topo_max, \
                    local_to_global_map,  \
                    lcolfrac, cell_id(I,J,nI,nJ))
    
            if lcolfrac:
                if count_col > 0:
                    global_field[J, I, :] = fsum/count_col
            else: 
                valid = count > 0
                global_field[J, I, valid] = fsum[valid]/count[valid]
            
    return global_field