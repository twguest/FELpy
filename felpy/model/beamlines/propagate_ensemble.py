#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "0.1.1"
__maintainer__ = "Trey Guest"
__email__ = "twguest@students.latrobe.edu.au"
__status__ = "Developement"
"""

from os import listdir
from felpy.utils.os_utils import mkdir_p

from felpy.model.beamlines.exfel_spb.methods import get_beamline_object
from felpy.model.core.wavefront import Wavefront
from felpy.model.tools import propagation_parameters, scale
from felpy.model.src.coherent import construct_SA1_pulse
from wpg.optical_elements import Drift
from wpg.wpg_uti_wf import plot_intensity_map as plot_wfr
from multiprocessing import Pool, cpu_count
from functools import partial
from wpg import srwlib
def load_wfr(indir):
    wfr = Wavefront()
    wfr.load_hdf5(indir)
    return wfr

def propagate_from_ensemble(wfr_file, ensemble_dir, bl, ekev, sdir = None, scale_input = None,
                            VERBOSE = True):
    """
    propagate a single wavefont from an ensemble of wavefronts (stored in hdf5 files) 
    down a desired beamline. This script can be used in conjunction with jobscheduler
    
    :param wfr_file: llocation of .hdf5 file
    :param scale: bool or list of scaling parameters [nx, ny, fov] for scaling 
    of input wavefield (to be updated)
    """
    if VERBOSE:
        print("loading wavefront")
        print("ensemble dir: {}".format(ensemble_dir))
        print("file: {}".format(wfr_file))
        
    wfr = load_wfr(ensemble_dir + "/" + wfr_file)
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
    
    if VERBOSE: 
        print("wavefront loaded")
        print(wfr.srw_info())
        
    if scale_input is not None and type(scale) == list:
            scale(wfr, scale[0], scale[1], scale[2])
    elif scale:
            wfr = scale(wfr)
            if VERBOSE:
                print("wavefront scaled")
    
    if VERBOSE:
        print("propagating wavefront")
        print(bl)
        
    bl.propagate(wfr)
    
    if VERBOSE:
        print("wavefront propagated")
    
    if sdir is not None:
        mkdir_p(sdir)
        wfr.store_hdf5(sdir + wfr_file)
        
def propagate_ensemble(ensemble, bl, sdir = None, scale_input = None, VERBOSE = True):
    """
    propagate an ensemble of wavefronts (stored in hdf5 files) in ensemble_dir
    down a desired beamline. 
    
    Useful in the case when there is a lot of ensembles
    
    :param ensemble: list of strings pointing to hdf5 wavefront files
    :param scale: bool or list of scaling parameters [nx, ny, fov] for scaling 
    of input wavefield (to be updated)
    """
    for file in ensemble:
        

        wfr = load_wfr(file)
            
        if scale_input is not None and type(scale) == list:
            scale(wfr, scale[0], scale[1], scale[2])
        elif scale:
            scale(wfr)
        
        bl.propagate(wfr)
        
        if sdir is not None:
            mkdir_p(sdir)
            wfr.store_hdf5(sdir + file)
            
def propagate_ensemble_mpi(ensemble, bl, sdir = None, scale_input = None,
                       processes = cpu_count()//2):
    """
    Useful in the case of many small pulses
    
    Note that mpi may by limited by available RAM.
    """
    
    mpi = Pool(processes = processes)
    func = partial(propagate_ensemble, bl = bl, sdir = sdir,
                   scale_input = scale_input)

    mpi.map(func, (file for file in ensemble))
    
