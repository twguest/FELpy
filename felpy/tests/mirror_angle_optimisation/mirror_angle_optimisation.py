# -*- coding: utf-8 -*-
"""
loop to optimise the angle of the HOM1 mirror.
"""
import sys
import shutil

import numpy as np
import multiprocessing as mp

from felpy.model.beamlines.structure import BeamlineModel
from felpy.model.src.coherent import construct_spb_pulse
from felpy.analysis.energy_statistics import get_pulse_energy
from felpy.utils.os_utils import mkdir_p
from felpy.utils.vis_utils import animate

from wpg.wpg_uti_wf import plot_intensity_map as plot_wfr
from functools import partial

FOCUS = 'nano' ### applies for NKB beamline analysis
PROCESSES = mp.cpu_count()//2 ### number of cpus for multiprocessing
MIRROR = "HOM2" ### mirror angle to be optimised
N = 50 ### number of data points

sdir = "./{}/".format(MIRROR)


def core(ang, mirror_name):
    """
    Core function for mirror optimisation tests. Can be generalised for all
    mirrors. In future, may add arg method= for defining which metric to 
    optimise for
    
    :param ang: mirror_angle
    
    :return data: set of optimisation parameters and output metrics. 
    """
    
    mkdir_p(sdir)
    
    spb = BeamlineModel(VERBOSE = False)
    spb.mirror_profiles(toggle = "on", aperture  = True, overwrite = False)
    spb.adjust_mirror("HOM1", 5.0, new_ang = ang)
    spb.buildElements(focus = FOCUS)
    spb.buildBeamline(focus = FOCUS)
    spb.cropBeamline(spb.params['HOM1']['next_drift'])
    bl = spb.get_beamline()
    
    wfr = construct_spb_pulse(512, 512, 2, 5.0, 0.25)
    
    bl.propagate(wfr)
    plot_wfr(wfr, save = sdir + "{:.4f}.png".format(ang*1e3))
    nph = get_pulse_energy(wfr)[1]
    print("Mirror Angle: {} mrad Complete".format(ang*1e3
                                                 ))
    return (ang, nph)
    
    
if __name__ == '__main__':
    
    MIRRORS = ['HOM1', 'HOM2', 'NHE', 'NVE']
    
    for MIRROR in MIRRORS:    
    
        spb = BeamlineModel(VERBOSE = False)
        pool = mp.Pool(PROCESSES)
        r = pool.map(partial(core, mirror_name = MIRROR),
                 np.linspace(spb.params[MIRROR]['ang_min'],
                             spb.params[MIRROR]['ang_max'],
                             N))
        
        
        animate(sdir, sdir, "{}_mirror_rotation".format(MIRROR))
        
        np.save(r, sdir + "{}_mirror_flux_data".format(MIRROR))
        
        