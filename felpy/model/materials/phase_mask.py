#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "1.0.0"
__maintainer__ = "Trey Guest"
__email__ = "twguest@students.latrobe.edu.au"
__status__ = "Developement"
"""

import numpy as np

from wpg.srwlib import SRWLOptD as Drift
from felpy.model.src.coherent import construct_SA1_wavefront
from felpy.model.beamline.structure import propagation_parameters
from wpg.beamline import Beamline
from wpg.srwlib import srwl_opt_setup_surf_height_2d as OPD
from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity
from matplotlib import pyplot as plt
from felpy.model.materials.material_utils import add_extent

def pltPhase(wfr):
    phase = wfr.get_phase()[:,:,0]
    plt.imshow(phase, cmap = 'hsv')
    plt.show()

def phaseMask(phaseshift, extent, wav, _ang = 0, outdir = None, maskName = None):
    """
    :param phaseshift: 2D array of desired phase-shifts
    :param extent: extent of phase-mask array in realspace
    :param wav: radiation wavelength 
    """
    
    height_error = (phaseshift*wav)/(2*np.pi)
    height_error = add_extent(height_error, extent)

    if outdir is not None:
        
        if maskName is not None:
            outdir = outdir + maskName + ".dat"
        else:
            outdir = outdir + "phasemask.dat"
        
        np.savetxt(outdir, height_error)
                   
    return OPD(height_error,
               _dim = 'x',
               _ang = np.pi/2,
               _refl = 1,
               _x = 0, _y = 0)

if __name__ == "__main__":
    

    wfr = construct_SA1_wavefront(1024,1024,3, 1)
    pltPhase(wfr)
    plotIntensity(wfr)
    phaseshift = np.random.rand(200,200)*2*np.pi

    opd = phaseMask(phaseshift, [5e-03, 5e-03], wfr.params.wavelength*100)
    bl = Beamline()
    
    bl.append(opd,propagation_parameters(1,1,1,1))
    bl.append(Drift(2), propagation_parameters(1,1,1,1))
    bl.propagate(wfr)
    pltPhase(wfr)
    plotIntensity(wfr)
    
    