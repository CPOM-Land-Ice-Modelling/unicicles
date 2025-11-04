#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 09:04:03 2025

@author: ggslc
"""
import numpy as np
import unittest
import matplotlib.pyplot as plt
from downscale import interp_to_local,interp_to_surface
from upscale import mean_to_global_mec
from test_scream_data import scream_surface_test_data

#generate default 'scream' test data 
X, Y, XX, YY, xXY, yXY, x, y, xx, yy, \
        Ixy, Jxy, q_sftc_g, q_topo_g, usrf, mask,topo_max_ec = scream_surface_test_data()
  
#interpolate global fields for each elevation class
q_sftc_l = interp_to_local(q_sftc_g, xXY, yXY, xx, yy)
q_topo_l = interp_to_local(q_topo_g, xXY, yXY, xx, yy)

#project onto surface
sftc = mask*interp_to_surface(q_sftc_l, q_topo_l, usrf)
       
#upscale the downscaled temperature  
p_sftc_g = mean_to_global_mec(sftc, usrf, topo_max_ec, \
                            (Ixy, Jxy), q_sftc_g.shape)

    
class Test_mean_to_global_mec(unittest.TestCase):
    
    def test_max(self):
        self.assertAlmostEqual(np.max(acab_g), 0.7321309239418563)

    def test_mim(self):
        self.assertAlmostEqual(np.min(acab_g), -0.2499999999999878)

    def test_sumsq(self):
        self.assertAlmostEqual(np.sum(acab_g**2),73.97217312122652)
        
unittest.main()


if True:

    
    NEC = 5
    fig, axs = plt.subplots(4,3,figsize=(8,10))    
    fig.suptitle('Downscaling', fontsize=16)
    jndx = 0
    
    def temp_color_map(ax, x, y, T):
        return ax.pcolormesh(x, y, T,  vmin=-10,vmax=10,cmap='bwr')
    
    for NEC in [0,2,4]:
        #pc = axs.flat[jndx].pcolormesh(X, Y, q_sftc_g[:,:,NEC], vmin=-1,vmax=1,cmap='bwr_r')
        #fig.colorbar(pc,ax=axs.flat[jndx])
        #jndx += 1
        
        pc = temp_color_map(axs.flat[jndx], xXY, yXY, q_sftc_g[:,:,NEC])
        #fig.colorbar(pc,ax=axs.flat[jndx])
        axs.flat[jndx].axhline(np.min(y),color='k',lw=0.5)
        axs.flat[jndx].axhline(np.max(y),color='k',lw=0.5)
        axs.flat[jndx].axvline(np.min(x),color='k',lw=0.5)
        axs.flat[jndx].axvline(np.max(x),color='k',lw=0.5)
        axs.flat[jndx].set_title(f'global/coarse, ec = {NEC}')
        axs.flat[jndx].set_xticks([])
        axs.flat[jndx].set_yticks([])
        jndx += 1
        
    for NEC in [0,2,4]:      
        pc = temp_color_map(axs.flat[jndx], x, y, q_sftc_l[:,:,NEC])
        #fig.colorbar(pc,ax=axs.flat[jndx])
        axs.flat[jndx].set_title(f'local/fine, ec = {NEC}')
        axs.flat[jndx].set_xticks([])
        axs.flat[jndx].set_yticks([])
        jndx += 1
    
    pc = temp_color_map(axs.flat[jndx], x, y, sftc)
    #fig.colorbar(pc,ax=axs.flat[jndx])
    #axs.flat[jndx].contour(x,y,usrf,[0,500,1000,1500],colors=['k'])
    axs.flat[jndx].set_xlim([-180,180])
    axs.flat[jndx].set_ylim([-180,180])
    axs.flat[jndx].set_title('local/fine, at surface ')
    axs.flat[jndx].set_xticks([])
    axs.flat[jndx].set_yticks([])
    jndx += 1
    
    axs.flat[jndx].plot(q_topo_g.flat[:]-20.0, q_sftc_g.flat,'.',label='g_sftc_q',ms=1)
    axs.flat[jndx].plot(q_topo_l.flat[::77]+20, q_sftc_l.flat[::77],'.',label='q_sftc_l',ms=1)
    axs.flat[jndx].plot(usrf.flat[::77],sftc.flat[::77],'.',label='sftc',color='k',ms=1)
    axs.flat[jndx].legend()
    jndx += 1
 
    fig, axs = plt.subplots(2, 3, figsize=(12,6))
    nec = p_sftc_g.shape[2]
    for ec in range(0, nec):
        ax = axs.flat[ec]
        pc = ax.pcolormesh(p_sftc_g[:, :, ec], vmin=-1, vmax=1, cmap='bwr_r')
        ax.set_title(f'{topo_max_ec[ec]} < s < {topo_max_ec[ec+1]}')
        #ax.contour(usrf_g[:, :, ec], [topo_max_ec[ec], topo_max_ec[ec+1]])
        fig.colorbar(pc)

    fig, axs = plt.subplots(1, 1)
    sftc_vsum = np.sum(p_sftc_g[:, :, :], axis=2)
    pc = axs.pcolormesh(sftc_vsum, vmin=-1, vmax=1, cmap='bwr_r')
    fig.colorbar(pc)

