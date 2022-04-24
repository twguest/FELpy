# -*- coding: utf-8 -*-
import numpy as np

from felpy.utils.vis_utils import scatter_plot
from felpy.model.beamlines.exfel_spb.exfel_spb import Instrument
from felpy.model.beamlines.exfel_spb.methods import get_beamline_object
from felpy.model.src.coherent import construct_SA1_wavefront
from tqdm import tqdm
from wpg.wpg_uti_wf import plot_intensity_map
from felpy.model.materials.load_refl import load_refl, get_refl 

import multiprocessing as mpi
from felpy.utils.np_utils import memory_map, readMap

from functools import partial



def no_mirror(ekev = 5.0):
    wfr= construct_SA1_wavefront(512, 512, ekev, 0.25)
        
    bl = get_beamline_object(ekev, apertures = False, surface = False,
                             crop = ["d1", "d1"])
    
    bl.propagate(wfr)
    ii = wfr.get_intensity().sum()
    plot_intensity_map(wfr)
    
    return ii

 
def no_aperture(angle = 1e-03, ekev = 5.0):
    
   
    wfr= construct_SA1_wavefront(512, 512, ekev, 0.25)
    
    bl = get_beamline_object(ekev, apertures = False, surface = True,
                             crop = ["d1", "HOM1"], theta_HOM = angle)
    
    bl.propagate(wfr)
    print(wfr.get_intensity().sum())
    return wfr.get_intensity().sum()

def aperture(angle = 1e-03, ekev = 5.0):
    
   
    wfr= construct_SA1_wavefront(512, 512, ekev, 0.25)
    
    bl = get_beamline_object(ekev, apertures = True, surface = True,
                             crop = ["d1", "HOM1"], theta_HOM = angle)
    
    bl.propagate(wfr)
     
    return wfr.get_intensity().sum()


if __name__ == '__main__':
    
    from labwork.about import dCache
    from felpy.utils.os_utils import mkdir_p
    
    try:
        SDIR =   "./mirror_reflectivity/"
        mkdir_p(SDIR)

    except(FileNotFoundError):
        SDIR = input("Save Directory: ")
        mkdir_p(SDIR)

    
    #ii = no_mirror()
    
    energies = [5.0, 7.0, 9.0, 11.0, 12.0]
    angles = np.linspace(0, 10e-03, 35)
    
    noap = memory_map(SDIR + "mirror_refl_no_aperture", shape = (len(energies), len(angles),2))
    ap = memory_map(SDIR + "mirror_refl_aperture", shape = (len(energies), len(angles),2))
    
    cpus = mpi.cpu_count()//2
    print("Distributing to {} cpus".format(cpus))
    
    for a in tqdm(range(len(energies))):
        
        pool = mpi.Pool(processes = cpus)
        
        #ii = no_mirror(energy)
        func_a = partial(no_aperture, ekev = energies[a])
      
        res = pool.map(func_a, angles)
        noap[a, :, :] = np.array([angles, res]).T
        
        func_b = partial(aperture, ekev = energies[a])
        res = pool.map(func_b, tqdm(angles))
        ap[a, :, :] = np.array([angles, res]).T
        
        
    #data = no_aperture()
    data = readMap(SDIR + "mirror_refl_no_aperture", shape = (len(energies), len(angles),2))