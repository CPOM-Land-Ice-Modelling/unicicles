#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 30 16:13:44 2025

@author: ggslc
"""

import numpy as np
from functools import partial


def subset_map(local_to_global_map, I, J):
    I_lg, J_lg = local_to_global_map
    return (I_lg == I) & (J_lg == J)


def subset_ec(topo, topo_max, ec):
    return (topo < topo_max[ec+1]) & (topo >= topo_max[ec])


def mean_to_global_cell_mec(field, topo, topo_max, local_to_global_map,
                            lcolfrac, I, J):

    indx = subset_map(local_to_global_map, I, J)

    nec = len(topo_max) - 1  # 0,z1,z2,,,,zn
    field_sum = np.zeros(nec)
    field_count = np.zeros(nec)
    ts = topo[indx]
    tf = field[indx]
    field_count_col = len(tf) + 1.0e-10
    for ec in range(0, nec):
        kndx = subset_ec(ts, topo_max, ec)
        t = tf[kndx]
        field_sum[ec] = np.sum(t)
        field_count[ec] = len(t) + 1.0e-10

    field_sum /= (field_count_col if lcolfrac else field_count)

    # conservation check
    # lsum = np.sum(field[indx])
    # gsum = np.sum(field_sum * field_count)

    return field_sum


def mean_to_global_mec(field, topo, topo_max, local_to_global_map,
                       global_shape, lcolfrac=False):

    nI, nJ, __ = global_shape
    global_field = np.zeros(global_shape)

    for I in range(0, nI):
        f = partial(mean_to_global_cell_mec, field, topo, topo_max,
                    local_to_global_map, lcolfrac, I)
        for J in range(0, nJ):
            global_field[I, J, :] = f(J)

    return global_field