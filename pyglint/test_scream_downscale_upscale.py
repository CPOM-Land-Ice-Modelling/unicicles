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
from transformation import cell_id

#generate default 'scream' test data 
X, Y, XX, YY, xXY, yXY, x, y, xx, yy, \
       local_to_global_map, sftc_g, topo_g, usrf, \
           mask,topo_max_ec = \
               scream_surface_test_data(local_grid_shape = (480,512),
                                        global_grid_shape = (18,16),
                                        global_min_elev = 100.0,
                                        global_max_elev = 1700, 
                                        global_grid_angle=np.pi/2,
                                        n_elev = 5)
  

#n = 0
#usrf = usrf*0 + 0.5*(topo_max_ec[n] + topo_max_ec[n+1]) 
#mask = mask*0 + 1
#sftc_g = sftc_g*0 + 5
  
#interpolate global fields for each elevation class
sftc_l = interp_to_local(sftc_g, xXY, yXY, xx, yy, order=0)
topo_l = interp_to_local(topo_g, xXY, yXY, xx, yy, order=0)

#project onto surface
sftc = mask*interp_to_surface(sftc_l, topo_l, usrf, lapse=0.0065)
       
#upscale the downscaled temperature  
p_sftc_g = mean_to_global_mec(sftc, usrf, topo_max_ec, \
                           local_to_global_map, sftc_g.shape)

    
class TestDownscale(unittest.TestCase):    
    
    def test_max(self):
        self.assertAlmostEqual(np.max(sftc), 13.44931497764644)
        
    def test_min(self):
        self.assertAlmostEqual(np.min(sftc), -9.47081230980631)
        
    def test_sumsq(self):
        self.assertAlmostEqual(np.sum(sftc**2),2933223.093176326)    
        
def test_sumsq(self):
        self.assertAlmostEqual(np.sum(p_sftc_g**2),10905.463676944699)
        
    
class TestUpscale(unittest.TestCase):
    
    def test_max(self):
        self.assertAlmostEqual(np.max(p_sftc_g), 11.829876400292084)

    def test_min(self):
        self.assertAlmostEqual(np.min(p_sftc_g), -8.906833236493016)

    def test_sumsq(self):
        self.assertAlmostEqual(np.sum(p_sftc_g**2),10905.463676944699)
        
        
unittest.main()


if True:

    
    fig, axs = plt.subplots(4,3,figsize=(8,10))    
    fig.suptitle('Downscaling', fontsize=16)
    jndx = 0
    
    def temp_color_map(ax, x, y, T):
        return ax.pcolormesh(x, y, T,  vmin=-10,vmax=10,cmap='bwr')
    
    for NEC in [0,1,2,3,4]:
        #pc = axs.flat[jndx].pcolormesh(X, Y, sftc_g[:,:,NEC], vmin=-1,vmax=1,cmap='bwr_r')
        #fig.colorbar(pc,ax=axs.flat[jndx])
        #jndx += 1
        
        pc = temp_color_map(axs.flat[jndx], xXY, yXY, sftc_g[:,:,NEC])
        #fig.colorbar(pc,ax=axs.flat[jndx])
        axs.flat[jndx].axhline(np.min(y),color='k',lw=0.5)
        axs.flat[jndx].axhline(np.max(y),color='k',lw=0.5)
        axs.flat[jndx].axvline(np.min(x),color='k',lw=0.5)
        axs.flat[jndx].axvline(np.max(x),color='k',lw=0.5)
        axs.flat[jndx].set_title(f'global/coarse, ec = {NEC}')
        axs.flat[jndx].set_xticks([])
        axs.flat[jndx].set_yticks([])
        jndx += 1
        
    for NEC in [0,1,2,3,4]:      
        pc = temp_color_map(axs.flat[jndx], x, y, sftc_l[:,:,NEC])
        #fig.colorbar(pc,ax=axs.flat[jndx])
        axs.flat[jndx].set_title(f'local/fine, ec = {NEC}')
        axs.flat[jndx].set_xticks([])
        axs.flat[jndx].set_yticks([])
        jndx += 1
    
    pc = temp_color_map(axs.flat[jndx], x, y, sftc)
    #fig.colorbar(pc,ax=axs.flat[jndx])
    #axs.flat[jndx].contour(x,y,usrf,[0,500,1000,1500],colors=['k'])
    axs.flat[jndx].set_title('local/fine, at surface ')
    axs.flat[jndx].set_xticks([])
    axs.flat[jndx].set_yticks([])
    jndx += 1
    
    axs.flat[jndx].plot(topo_g.flat[:]-20.0, sftc_g.flat,'.',label='sftc_g',ms=1)
    axs.flat[jndx].plot(topo_l.flat[::77]+20, sftc_l.flat[::77],'.',label='sftc_l',ms=1)
    axs.flat[jndx].plot(usrf.flat[::77],sftc.flat[::77],'.',label='sftc',color='k',ms=1)
    axs.flat[jndx].legend()
    jndx += 1
 
    
# %% Upscaling figure
    fig, axs = plt.subplots(4, 3, figsize=(8,10))
    fig.suptitle('Upscaling', fontsize=16)
    nec = p_sftc_g.shape[2]
    jndx = 0
    
    
    M, N = sftc_g.shape[0:2]
    J = np.arange(0,M)
    I = np.arange(0,N)
    II, JJ = np.meshgrid(I, J)
    CC = cell_id(II, JJ, N, M)
    
    replace_where = np.where(np.isfinite(p_sftc_g),p_sftc_g, sftc_g)
    err = replace_where-sftc_g
    
    
    
    for ec in range(0, nec):
        
        pc = temp_color_map(axs.flat[jndx], xXY , yXY,  p_sftc_g[:, :, ec])
        axs.flat[jndx].set_title(f'{topo_max_ec[ec]} < s < {topo_max_ec[ec+1]}')
        axs.flat[jndx].set_xticks([])
        axs.flat[jndx].set_yticks([])
        axs.flat[jndx].axhline(np.min(y),color='k',lw=0.5)
        axs.flat[jndx].axhline(np.max(y),color='k',lw=0.5)
        axs.flat[jndx].axvline(np.min(x),color='k',lw=0.5)
        axs.flat[jndx].axvline(np.max(x),color='k',lw=0.5) 
        axs.flat[jndx].contour(xXY , yXY, II,I+0.5,linewidths=0.5)
        axs.flat[jndx].contour(xXY , yXY, JJ,J+0.5,linewidths=0.5)
        axs.flat[jndx].set_xlim(-300,300)
        axs.flat[jndx].set_ylim(-300,300)
       
        
        
        jndx += 1
        
      
        pc = temp_color_map(axs.flat[jndx], xXY , yXY, err[:,:,ec])
        axs.flat[jndx].contour(xXY , yXY, II,I+0.5,linewidths=0.5)
        axs.flat[jndx].contour(xXY , yXY, JJ,J+0.5,linewidths=0.5)
        axs.flat[jndx].set_xlim(-300,300)
        axs.flat[jndx].set_ylim(-300,300)
        #axs.flat[jndx].set_ylim(-250,250)
        
        axs.flat[jndx].set_title(f'{topo_max_ec[ec]} < s < {topo_max_ec[ec+1]}')
        axs.flat[jndx].set_xticks([])
        axs.flat[jndx].set_yticks([])
        axs.flat[jndx].axhline(np.min(y),color='k',lw=0.5)
        axs.flat[jndx].axhline(np.max(y),color='k',lw=0.5)
        axs.flat[jndx].axvline(np.min(x),color='k',lw=0.5)
        axs.flat[jndx].axvline(np.max(x),color='k',lw=0.5)
        
        
        #ax.contour(usrf_g[:, :, ec], [topo_max_ec[ec], topo_max_ec[ec+1]])
        #fig.colorbar(pc)
        jndx += 1
    
    #sftc_vsum = np.sum(p_sftc_g[:, :, :], axis=2)
    #pc = axs.flat[jndx].pcolormesh(sftc_vsum, vmin=-1, vmax=1, cmap='bwr_r')
    #fig.colorbar(pc)

