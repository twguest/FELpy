#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "1.0.1"
__maintainer__ = "Trey Guest"
__email__ = "twguest@students.latrobe.edu.au"
__status__ = "Developement"
"""

import numpy as np

from wpg.srwlib import SRWLOptD as Drift
from felpy.model.source.coherent import construct_SA1_wavefront
from felpy.model.core.beamline import Beamline
from wpg.srwlib import srwl_opt_setup_surf_height_2d as OPD
from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity
from matplotlib import pyplot as plt
from felpy.model.materials.material_utils import add_extent

def plot_phase(wfr):
    phase = wfr.get_phase()[:,:,0]
    plt.imshow(phase, cmap = 'hsv')
    plt.show()



def phase_mask(phase_shift, extent, wav, _ang = 0, outdir = None, maskName = None):
    """
    :param phase_shift: 2D array of desired phase-shifts
    :param extent: extent of phase-mask array in realspace
    :param wav: radiation wavelength 
    """
    
    height_error = (phase_shift*wav)/(2*np.pi)
    height_error = add_extent(height_error, extent)

    if outdir is not None:
        
        if maskName is not None:
            outdir = outdir + maskName + ".dat"
        else:
            outdir = outdir + "phase_mask.dat"
        
        np.savetxt(outdir, height_error)
                   
    return OPD(height_error,
               _dim = 'x',
               _ang = 0,
               _refl = 1,
               _x = 0, _y = 0)

    
    