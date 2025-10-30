#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 09:04:03 2025

@author: ggslc
"""
import numpy as np
import unittest
import matplotlib.pyplot as plt
from mean_to_global_mec import mean_to_global_mec

# test data / local grid. centers.
dx = 1.0
m, n = 1024, 1024

x = np.linspace(-m/2+dx/2, m/2-dx/2, m)
y = np.linspace(-n/2+dx/2, n/2-dx/2, n)

xx, yy = np.meshgrid(x, y)

# test surface
usrf = 0.7*np.cos(14*xx/m) + 0.3*np.sin(12*yy/n)
usrf = np.where(xx**2+yy**2 < (m/3)**2, 2000.0*np.abs(usrf), 0)

ELA = 500
acab = 0.5 * (usrf - ELA) / 1000

M, N = 32, 32 
L = np.max(x) - np.min(x)
W = np.max(y) - np.min(y)
I = np.int32( M/2* (xx-np.min(x)+yy-np.min(y))/L )
J = np.int32( N/2 + N/2* (yy-np.min(y)-(xx-np.min(x)))/L )

# I = np.int32(M * (xx-np.min(x)) / L)
# J = np.int32(N * (yy-np.min(y)) / W)

# test data, elevation classes
nec = 6
topo_max_ec = np.linspace(0, 2400, nec+1)

acab_g = mean_to_global_mec(acab, usrf, topo_max_ec, (J, I), (N, M, nec))


class Test_mean_to_global_mec(unittest.TestCase):
    

    def test_max(self):
        self.assertAlmostEqual(np.max(acab_g), 0.7321309239418563)

    def test_mim(self):
        self.assertAlmostEqual(np.min(acab_g), -0.2499999999999878)

    def test_sumsq(self):
        self.assertAlmostEqual(np.sum(acab_g**2),73.97217312122652)
        
unittest.main()


if False:

    fig, axs = plt.subplots(2, 3, figsize=(12,6))

    for ec in range(0, nec):
        ax = axs.flat[ec]
        pc = ax.pcolormesh(acab_g[:, :, ec], vmin=-1, vmax=1, cmap='bwr_r')
        ax.set_title(f'{topo_max_ec[ec]} < s < {topo_max_ec[ec+1]}')
        #ax.contour(usrf_g[:, :, ec], [topo_max_ec[ec], topo_max_ec[ec+1]])
        fig.colorbar(pc)

    fig, axs = plt.subplots(1, 1)
    acab_vsum = np.sum(acab_g[:, :, :], axis=2)
    pc = axs.pcolormesh(acab_vsum, vmin=-1, vmax=1, cmap='bwr_r')
    fig.colorbar(pc)

    fig, axs = plt.subplots(1, 1)
    pc = axs.pcolormesh(acab, vmin=-1, vmax=1, cmap='bwr_r')
    fig.colorbar(pc)

    fig, axs = plt.subplots(1, 1)
    pc = axs.pcolormesh(usrf, vmin=0, vmax=2500, cmap='jet')
    fig.colorbar(pc)