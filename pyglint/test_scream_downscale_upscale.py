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
               scream_surface_test_data(local_grid_shape = (512,511),
                                        global_grid_shape = (31,33),
                                        global_min_elev = 50.0,
                                        global_max_elev = 1700, 
                                        global_grid_angle=np.pi/6,
                                        n_elev = 5)
  

n = 0
#usrf = usrf*0 + 0.5*(topo_max_ec[n] + topo_max_ec[n+1]) 
#mask = mask*0 + 1
#sftc_g = sftc_g*0 + 5
  
#interpolate global fields for each elevation class
sftc_l = interp_to_local(sftc_g, xXY, yXY, xx, yy, order=1)
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


def temp_color_map(ax, xc, yc, T,cmap='plasma'):
    
    ax.set_aspect('equal')
    pc = ax.pcolormesh(xc, yc, T,  vmin=-10,vmax=10,cmap=cmap)

    return pc

 
    
# %% Downscaling figure
if True:
    
    fig = plt.figure(figsize=(8,11))    
    fig.suptitle('Downscaling', fontsize=14)
    
    margin = 0.025
    aspect = fig.get_figwidth()/fig.get_figheight()
    step = 2
    nec = p_sftc_g.shape[2]
    nec_range = range(0, nec, step)
    width = (1.0 - (1+len(nec_range))*margin)/len(nec_range)
    height = width*aspect
    for i,ec in enumerate(range(0, nec, step)):
        ax = fig.add_axes(((1+i)*margin + i*width, 1 - (3*margin+height), width,height))
        pc = temp_color_map(ax, xXY, yXY, sftc_g[:,:,ec])
        _ = [ax.axhline(yl,color='k',lw=0.5) for yl in [np.min(y),np.max(y)]]
        _ = [ax.axvline(xl,color='k',lw=0.5) for xl in [np.min(x),np.max(x)]]
        ax.set_title(f'global, {topo_max_ec[ec]} < s < {topo_max_ec[ec+1]}')
        ax.set_xticks([])
        if (ec > 0 ):
            ax.set_yticks([])
        
        ax = fig.add_axes(((1+i)*margin + i*width, 1 - (5*margin+2*height), width,height))
        pc = temp_color_map(ax, x, y, sftc_l[:,:,ec])
        ax.set_title(f'local, {topo_max_ec[ec]} < s < {topo_max_ec[ec+1]}')
        if (ec > 0 ):
            ax.set_yticks([])
            
    #large, tsfc on surface
    ax = fig.add_axes((margin, margin, 1.75*width,1.75*height))
    pc = temp_color_map(ax, x, y, sftc)
    fig.colorbar(pc,ax=ax,shrink = 0.5, pad=margin)
    ax.set_title('local/fine, at surface ')   
    
    #large, distribution
    ax = fig.add_axes((4*margin + 1.75*width, 2.25*margin, 2*margin+width,1.4*height))
    ax.plot(sftc_g.flat,topo_g.flat[:]-20.0, 's',label='sftc_g (z+20m)',ms=2)
    ax.plot(sftc_l.flat[::77], topo_l.flat[::77]+20,'^',label='sftc_l (z-20m)',ms=2)
    ax.plot(sftc.flat[::77],usrf.flat[::77],'o',label='sftc',color='k',ms=2)
    ax.set_yticks(topo_max_ec,topo_max_ec/1000 )
    ax.grid()
    ax.set_xlabel(r'$T$ ($^\circ$C)')
    ax.set_ylabel(r'$z$ (km)')
    ax.legend(loc='upper right')
    
            
# %% Upscaling figure
if True:
  
    fig, axs = plt.subplots(3, 3, figsize=(8,11)) 
    fig.suptitle('Upscaling', fontsize=14)
     
    replace_where = np.where(np.isfinite(p_sftc_g),p_sftc_g, sftc_g)
    err = np.where(np.isfinite(p_sftc_g),p_sftc_g-sftc_g, 0.0)
    for i,ec in enumerate(range(0, nec, step)):
        for j, ffn in  enumerate(zip([sftc_g,p_sftc_g,err],['global','upscaled','error'])):
           ax = axs[j,i]
           f, fn = ffn
           pc = temp_color_map(ax, xXY, yXY, f[:,:,ec],cmap='bwr')
           _ = [ax.axhline(yl,color='k',lw=0.5) for yl in [np.min(y),np.max(y)]]
           _ = [ax.axvline(xl,color='k',lw=0.5) for xl in [np.min(x),np.max(x)]]
           ax.set_title(f'{fn}, ec = {ec}')
           ax.set_xticks([])
           fig.colorbar(pc,ax=ax,orientation='horizontal')
           if (ec > 0 ):
               ax.set_yticks([])
           
         




