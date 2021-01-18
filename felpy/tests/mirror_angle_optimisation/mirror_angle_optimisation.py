# -*- coding: utf-8 -*-
"""
loop to optimise the angle of the HOM1 mirror.
"""
import sys
import numpy as np
import multiprocessing as mp

from felpy.model.beamlines.structure import BeamlineModel
from felpy.model.src.coherent import construct_pulse
from felpy.analysis.energy_statistics import get_pulse_energy
from felpy.utils.os_utils import mkdir_p

from wpg.wpg_uti_wf import plot_intensity_map as plot_wfr
from functools import partial

FOCUS = 'nano' ### applies for NKB beamline analysis
PROCESSES = mp.cpu_count()//2 ### number of cpus for multiprocessing
MIRROR = "HOM1" ### mirror angle to be optimised
N = 2 ### number of data points

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
    spb.mirror_profiles(toggle = "on", aperture = True, overwrite = False)
    spb.adjust_mirror("HOM1", 5.0, new_ang = ang)
    print(1)
    spb.buildElements(focus = FOCUS)
    spb.buildBeamline(focus = FOCUS)
    spb.cropBeamline(spb.params['HOM1']['next_drift'])
    print(2)
    bl = spb.get_beamline()
    
    wfr = construct_pulse(512,512,2)
    
    bl.propagate(wfr)
    plot_wfr(wfr, save = sdir + "{:.6f}.png".format(ang*1e3))
    #nph = get_pulse_energy(wfr)[1]
    print("Mirror Angle: {} mrad Complete".format(ang*1e3
                                                  ))
    
    
    
if __name__ == '__main__':



    spb = BeamlineModel(VERBOSE = False)
    pool = mp.Pool(PROCESSES)
    pool.map(partial(core, mirror_name = MIRROR),
             np.linspace(spb.params[MIRROR]['ang_min'],
                         spb.params[MIRROR]['ang_max'],
                         N))
    