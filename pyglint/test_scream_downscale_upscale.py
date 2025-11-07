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
from test_scream_data import scream_surface_test_data, LAPSE_RATE
from transformation import cell_id, missing



def temp_color_map(ax, xc, yc, T,cmap='bwr'):
    
    ax.set_aspect('equal')
    pc = ax.pcolormesh(xc, yc, T,  vmin=-10,vmax=10,cmap=cmap)

    return pc

def downscale_fig(sftc_g, sftc_l, sftc, error,
                  topo_g, topo_l, usrf, 
                  xXY, yXY, x, y,  topo_max_ec):
    
    fig = plt.figure(figsize=(8,11))    
    fig.suptitle('Downscaling', fontsize=14)
    
    margin = 0.025
    aspect = fig.get_figwidth()/fig.get_figheight()
    step = 2
    nec = sftc_g.shape[2]
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
    ax.plot(sftc.flat[::77],usrf.flat[::77],'o',label='sftc (z)',color='k',ms=2)
    ax.plot(10*error.flat[::77],usrf.flat[::77],'o',label='10*error (z)',color='r',ms=2)
    ax.set_yticks(topo_max_ec,topo_max_ec/1000 )
    ax.grid()
    ax.set_xlabel(r'$T$ ($^\circ$C)')
    ax.set_ylabel(r'$z$ (km)')
    ax.legend(loc='upper right')
 
    
 
    
def upscale_fig(sftc_g,up_sftc_g,upscale_error,xXY, yXY, x, y):    
    
    fig, axs = plt.subplots(3, 3, figsize=(8,11)) 
    fig.suptitle('Upscaling', fontsize=14)
     
    nec = sftc_g.shape[2]
    step = 2
    for i,ec in enumerate(range(0, nec, step)):
        for j, ffn in  enumerate(zip([sftc_g,up_sftc_g,upscale_error],
                                     ['global','upscaled','error'])):
           ax = axs[j,i]
           f, fn = ffn
           pc = temp_color_map(ax, xXY, yXY, f[:,:,ec])
           _ = [ax.axhline(yl,color='k',lw=0.5) for yl in [np.min(y),np.max(y)]]
           _ = [ax.axvline(xl,color='k',lw=0.5) for xl in [np.min(x),np.max(x)]]
           ax.set_title(f'{fn}, ec = {ec}')
           ax.set_xticks([])
           fig.colorbar(pc,ax=ax,orientation='horizontal')
           if (ec > 0 ):
               ax.set_yticks([])
           


def scream_test(*args, **kwargs):
    

    #generate  'scream' test data 
    X, Y, XX, YY, xXY, yXY, x, y, xx, yy, \
           local_to_global_map, sftc_g, topo_g, usrf, sftc_known, \
               mask,topo_max_ec = scream_surface_test_data(*args, **kwargs)
 
    interp_lapse_rate = kwargs.get("interp_lapse_rate",LAPSE_RATE)
    
    #interpolate global fields for each elevation class
    sftc_l = interp_to_local(sftc_g, xXY, yXY, xx, yy, order=1)
    topo_l = interp_to_local(topo_g, xXY, yXY, xx, yy, order=0)
    
    #project onto surface
    sftc = np.where(mask, interp_to_surface(sftc_l, topo_l, usrf, 
                                            lapse = interp_lapse_rate), np.nan)
        
    #compare stfc to sftc_known
    downscale_error =  np.where(mask, (sftc - sftc_known), 0.0)  
    
    #upscale the downscaled temperature  
    up_sftc_g = mean_to_global_mec(sftc, usrf, mask, topo_max_ec, \
                               local_to_global_map, sftc_g.shape)
    
        
        
    upscale_error = np.where(np.isfinite(up_sftc_g),up_sftc_g-sftc_g, 0.0)
    
    if plot:
        downscale_fig(sftc_g, sftc_l, sftc,  downscale_error,
                      topo_g, topo_l, usrf, 
                      xXY, yXY, x, y, topo_max_ec)
    
        upscale_fig(sftc_g, up_sftc_g, upscale_error, 
                xXY, yXY, x, y)
    
    
    return downscale_error, upscale_error
    

# running the test codes here rather than in TestScream
# because the profiler works better that way...

plot = False

# case 0 : defaults
down_err_0, up_err_0 = scream_test() 

# case 0, flat surface, atmosphere uniform in (x,y)
# z = 150 is the middle of the lowest elevation class
down_err_1, up_err_1 = scream_test(z_freeze_var=0.0, 
                                   usrf_max=0.0, usrf_min=150.0) 

    
class TestScream(unittest.TestCase):    
 
    def test_down_0_mean(self):
        self.assertLess(np.abs(np.mean(down_err_0)),1.0e-3)
        
    def test_down_0_std(self):
        self.assertLess(np.std(down_err_0),0.03)
 
    def test_up_0_mean(self):
        self.assertLess(np.abs(np.mean(up_err_0)),1.0e-2)
    
    def test_up_0_std(self):
        self.assertLess(np.std(down_err_0),0.3)
 
    
    def test_down_1(self):
        # case 1 downscale should be exact  
        self.assertLess(np.max(np.abs(down_err_1)), 1.0e-10)
        
    def test_up_1(self):   
        #case 1 upscale should be exact
        self.assertLess(np.max(np.abs(up_err_1)), 1.0e-10)   
    
#class TestUpscale(unittest.TestCase):
    


    
unittest.main()



 
    
  

         




