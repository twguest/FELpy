#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 11:06:56 2020

@author: twguest

"""


###############################################################################
import sys
sys.path.append("/opt/WPG/") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/WPG") # DESY MAXWELL PATH

sys.path.append("/opt/spb_model") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/spb_model") # DESY MAXWELL PATH
###############################################################################
###############################################################################

import numpy as np
from copy import deepcopy

from felpy.model.src.coherent import construct_SA1_pulse
from felpy.model.instrument import Instrument

from matplotlib import pyplot as plt
from felpy.model.materials.load_refl import load_refl, get_refl

from tqdm import tqdm
from wpg.wpg_uti_wf import get_profile_1d
from felpy.utils.os_utils import mkdir_p
from wpg.wpg_uti_wf import get_axis

def beamline_setup(mode, ekev, ang = 3.3e-03):
    """
    Setup beamline, and adjust mirror parameters.
    Constructs two beamlines relating to the drift to the entry plane, and
    through the exit plane of HOM1
    
    :param mode: 'pre' or 'post' returns a beamline for propagation to the 
                entry plane of the optic (pre) and from the entry plane 
                to the exit plane (post)
    :param ekev: radiation energy [keV]
    :param ang: fixed mirror angle (defines mirror aperture)
    """
    spb = Instrument(VERBOSE = False)
 
    spb.adjust_mirror("HOM1", ekev, ang)
    spb.build_elements(focus = "nano")
<<<<<<< HEAD:felpy/tests/aperture_test/aperture_test.py
    spb.build_beamline(focus = "nano")
=======
    spb.buildBeamline(focus = "nano")
>>>>>>> 108cfb9b6fc97d3841ee1db54862523eee5b184e:felpy/examples/aperture_test/aperture_test.py
    
    if mode == 'pre': ### OPTION FOR PROPAGATING TO THE MIRROR SURFACE
        spb.crop_beamline(element1 = "d1") ### ie we propagate to the plane prior
    elif mode == 'post':
        spb.crop_beamline(element1 = "HOM1", element2 = "HOM1") ### ie we propagate to the plane prior
        
    bl = spb.get_beamline()
    
    return bl

def line_plot(pre, post, axis, sdir = None):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
           
    
    ax.plot(axis*1e3, pre, '--')
    
    for ix in post:
        ax.plot(axis*1e3, ix)
   
    
    ax.set_xlabel("Position (mm)")
    ax.set_ylabel("Intensity (W/$mm^2$)")

    
    plt.legend(["Wavefield in Entry Plane", "Wavefield in Exit-Plane"], loc = 1)
    plt.show()
    if sdir is not None:
        fig.savefig(sdir)
    
    


def get_unobstructed_profile(ekev, q):
    """
    Return one-dimensional profiles of the source in the entry-plane of the aperture
    
    :param ekev: radiation energy [keV]
    :param q: electron beam charge [nC]
    
    :returns ix: horizontal intensity profile
    :returns iy: vertical intensity profile
    :returns ax: spatial axis along x profilefrom tqdm import tqdm
    :returns ay: spatial axis along y profile
    """

    
    wfr = construct_SA1_pulse(1024, 1024, 2, ekev, q)
    bl = beamline_setup(mode = 'pre', ekev = ekev)

    bl.propagate(wfr)    
    ix, iy = get_profile_1d(wfr)

    return ix, iy, get_axis(wfr, axis = 'x'), get_axis(wfr, axis = 'y')


def get_obstructed_profile(wfr, bl):
    """
    Return one-dimensional profiles of the source in the exit-plane of the aperture
    
    :param ekev: radiation energy [keV]
    :param q: electron beam charge [nC]
    
    :returns ix: horizontal intensity profile
    :returns iy: vertical intensity profile
    """
    
  
    bl.propagate(wfr)    
    ix,iy = get_profile_1d(wfr)

    return ix, iy


def mirror_acceptance_test(ekev, q, A):
    
    iwfr = construct_SA1_pulse(1024, 1024, 2, ekev, q)
    
    ix, iy, axis_x, axis_y = [get_obstructed_profile(deepcopy(iwfr), beamline_setup(mode = 'post',
                                                      ekev = ekev,
                                                      ang = ang)) for ang in tqdm(A)]
    ux, uy = get_unobstructed_profile(ekev, q)
    line_plot(ix, ux, axis_x)
    

    
if __name__ == '__main__':
    
    ix, iy = mirror_acceptance_test(ekev = 5.0, q = 0.25, A = np.linspace(1.1, 3.6, 2))
