#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 16:10:46 2020

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

import time
from copy import copy
import numpy as np
from matplotlib import pyplot as plt
from model.src.coherent import coherentSource
from model.beamline.structure import buildBeamline, cropBeamline, load_params, adjustHOM1
from model.materials.load_refl import load_refl, get_refl



def define_wfr(ekev):
    """
    defines the wavefront in the plane prior to the mirror ie., after d1
    
    :param ekev: energy of the source
    """
    params = load_params()
    bl = buildBeamline(params)
    cropBeamline(bl, element1 = "drift1", element2 = "drift1")
    wfr = coherentSource(1048, 1048, ekev, 1)
    bl.propagate(wfr)
    
    return wfr

def getTransmission(ekev, ang):

    params = load_params()
    
    refl_data = load_refl()
    refl, ang = get_refl(refl_data, ekev, ang = ang, limits = [1.1e-03, 3.6e-03])
    
    print("reflectivity: {}".format(refl))
    print("incidence angle: {} mrad".format(ang*1e3))
    
    adjustHOM1(params, refl, ang)

    bl = buildBeamline(params)
    adjustHOM1(params, refl, ang)
    
    cropBeamline(bl, element1 = "HOM1", element2 = "HOM1")
    
    wfr = define_wfr(ekev)
    
    i = np.sum(wfr.get_intensity())
    bl.propagate(wfr)
    f = np.sum(wfr.get_intensity())

    transmission = f/i
    
    return transmission, refl, ang

def plotTransAng(ekev, dat_trans, dat_ang):
    
    refl = load_refl()
    
    angs = copy(refl[1:,0])
    
    e_idx = (np.abs(refl[0,1:] -ekev*1e3)).argmin() + 1
    refl = refl[:, e_idx][1: ]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.plot(angs, refl, '--')
    ax.scatter(dat_ang, dat_trans, c = 'r', marker = 'x')
    plt.yscale('log')
    
    ax.set_title("HOM1 Mirror Reflectivity")
    ax.set_xlabel("Mirror Grazing Angle (mrad)")
    ax.set_ylabel("Transmitted Intensity ($I_{f}/I_{i}$)")
    
    plt.legend(["Tabulated Data","Transmitted Intensity"])
    plt.show()
    
    estr = str(ekev).replace(".","-")
    fig.savefig("refl_test/mirror_refl_{}kev".format(estr))
    
def testTransmissionAngles(ekev, n = 100):
    
 
    ang_range = np.linspace(1e-03,15e-03, n)
    
    dat_trans = []
    dat_refl = []
    dat_ang = []
    
    for ang in ang_range:
        trans, refl, ang = getTransmission(ekev, ang)    
        
        dat_trans.append(trans)
        dat_refl.append(refl)
        dat_ang.append(ang)
        
    estr = str(ekev).replace(".","-")
    np.save("refl_test/data_trans_{}kev.npy".format(estr), dat_trans)    
    np.save("refl_test/data_ang_{}kev.npy".format(estr), dat_refl)   
    
    plotTransAng(ekev, dat_trans, dat_ang)

if __name__ == '__main__':
    
    energies = [6.0, 9.2, 12.0]
    
    for ekev in energies:
        testTransmissionAngles(ekev, n = 25)