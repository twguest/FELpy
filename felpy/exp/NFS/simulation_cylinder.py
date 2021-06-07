# -*- coding: utf-8 -*-
 
from felpy.model.materials.phase_mask import phase_mask
import numpy as np
 
from felpy.utils.opt_utils import ekev2wav
from felpy.exp.NFS.cylinder_phase_mask import phase 
from felpy.model.src.coherent import construct_gaussian

from wpg.wpg_uti_wf import plot_intensity_map as plot_ii
from wpg.optical_elements import Drift
from felpy.model.core.beamline import Beamline 
from felpy.model.tools import propagation_parameters

r1 = 10e-02
r2 = 5e-02

delta = 1e-04
ekev = 25
wav = ekev2wav(ekev)

k = (np.pi*2)/wav

x = np.linspace(-r1/2, r1/2, 100)

y = np.ones([100,100])

phase_shift = phase(x, r1, r2, k, delta)
mask = y*phase_shift





exp_setup = {}
exp_setup['z1'] = [1,2]
exp_setup['z2'] = [3,2]


for itr in range(len(exp_setup['z1'])):
    bl = Beamline()
    
    src = construct_gaussian(512, 512, ekev, extent = [-10e-03, 10e-03, -10e-03, 10e-03],
                         sigX = 1e-03, sigY = 1e-03, divergence = np.arctan(13.5e-03/8))
    
    bl.append(Drift(exp_setup['z1'][itr]), propagation_parameters(1, 1, 1, 1, mode = 'quadratic'))
    bl.propagate(src)
    plot_ii(src)
    
    r = (3*src.get_fwhm()[0])
    dx = 250
    x = np.linspace(-r/2, r/2, dx)

    y = np.ones([dx,dx])
    
    phase_shift = phase(x, r1, r2, k, delta)
    mask = np.ones([dx,dx])*phase_shift
    opt = phase_mask(mask, [-r/2, r/2, -1, 1], wav)
        
    bl = Beamline()
    bl.append(opt, propagation_parameters(1,1,1,1,mode = 'fresnel'))
    bl.append(Drift(exp_setup['z2'][itr]), propagation_parameters(1, 1, 1, 1, mode = 'quadratic'))
    bl.propagate(src)
    
    plot_ii(src)