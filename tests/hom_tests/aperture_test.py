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
from copy import copy

from model.src.coherent import coherentSource
from model.beamline.structure import BeamlineModel

from matplotlib import pyplot as plt
from model.materials.load_refl import load_refl, get_refl

from tqdm import tqdm

def beamlineSetup(mode, ekev):
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
    spb = BeamlineModel(overwrite_mirrors = True)
    
    refl_data = load_refl()
    refl, ang = get_refl(refl_data, ekev, ang = 2.2e-03, limits = [1.1e-03, 3.6e-03])
    
    spb.adjustHOMs(refl, ang) ### SETTING _refl to 1 for this test
    spb.buildElements(focus = "micron")
    spb.buildBeamline(focus = "micron")
    
    if mode == 'pre': ### OPTION FOR PROPAGATING TO THE MIRROR SURFACE
        spb.cropBeamline(element1 = "drift1") ### ie we propagate to the plane prior
    elif mode == 'post':
        spb.cropBeamline(element1 = "HOM1", element2 = "HOM1") ### ie we propagate to the plane prior
        
    bl = spb.get_beamline()
    
    return bl

def line_plot(pre, post, xx, ekev, q, label = 'Horizontal', outdir = None):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
           
    ax.plot(xx*1e3, pre, '--')
    ax.plot(xx*1e3, post)
   
    
    ax.set_xlabel("Position (mm)")
    ax.set_ylabel("Intensity (W/$mm^2$)")
    ax.set_title("{} Mirror Acceptance at {} keV, {} nC".format(label, ekev, q))
    
    
    
    plt.legend(["Wavefield in Entry Plane", "Wavefield in Exit-Plane"], loc = 1)
    plt.show()
    if outdir is not None:
        fig.savefig(outdir)
    
    
   
def get_pre(ekev, q):
    """
    Return one-dimensional profiles of the source in the entry-plane of the aperture
    
    :param ekev: radiation energy [keV]
    :param q: electron beam charge [nC]
    
    :returns ix: horizontal intensity profile
    :returns iy: vertical intensity profile
    :returns ax: spatial axis along x profilefrom tqdm import tqdm
    :returns ay: spatial axis along y profile
    """
    wfr = coherentSource(1024, 1024, ekev, q)
    bl = beamlineSetup(mode = 'pre', ekev = ekev)
    bl.propagate(wfr)    
    ix,iy = copy(wfr.get_profile_1d())
    
    
    ax = np.linspace(wfr.params.Mesh.xMin, wfr.params.Mesh.xMax, wfr.params.Mesh.nx)
    ay = np.linspace(wfr.params.Mesh.yMin, wfr.params.Mesh.yMax, wfr.params.Mesh.ny)
    
    return wfr, ix, iy, ax, ay


def get_post(wfr, ekev, q):
    """
    Return one-dimensional profiles of the source in the exit-plane of the aperture
    
    :param ekev: radiation energy [keV]
    :param q: electron beam charge [nC]
    
    :returns ix: horizontal intensity profile
    :returns iy: vertical intensity profile
    """
    
    bl = beamlineSetup(mode = 'post', ekev = ekev)
    bl.propagate(wfr)    
    ix,iy = copy(wfr.get_profile_1d())

    return ix, iy

def testAcceptance(ekev, q, outdir = None):
    """
    Test mirror acceptance
    """
    wfr, pre_x, pre_y, ax, ay = get_pre(ekev, q)
    post_x, post_y = get_post(wfr, ekev, q)

    
    estr = str(ekev).replace(".", "-")
    qstr = str(q).replace(".", "-")
        
    line_plot(pre_x, post_x, ax, ekev, q, label = "Horizontal", outdir = outdir + "horizontal_{}keV_{}nC".format(estr, qstr))
    line_plot(pre_y, post_y, ay, ekev, q, label = "Vertical", outdir = outdir + "vertical_{}keV_{}nC".format(estr, qstr))
    
    np.save("aperture_test/data/preX_{}keV_{}nC.npy".format(estr,qstr), pre_x)
    np.save("aperture_test/data/preY_{}keV_{}nC.npy".format(estr,qstr), pre_y)
    np.save("aperture_test/data/postX_{}keV_{}nC.npy".format(estr,qstr), post_x)
    np.save("aperture_test/data/postY_{}keV_{}nC.npy".format(estr,qstr), post_y)
    
    
if __name__ == '__main__':
    from tqdm import tqdm
    
    ang = 2.2e-03 # fixed mirror angle [rad]
    ### NOTE THAT MIRROR ANGLE IS HARDCODED IN beamlineSetup fn.
    
    print("Testing Mirror Acceptance for Fixed Mirror Angle: {} mrad".format(ang))
    
    energies = [6.0, 7.0, 8.0, 9.0, 9.2, 10.0, 11.0, 12.0]
    charges = [0.1, 0.25, 0.5]
    
    for ekev in tqdm(energies):
        for q in charges:
            testAcceptance(ekev, q, outdir = "aperture_test/mirror_acceptance_")
    

    
    