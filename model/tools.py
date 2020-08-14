
###############################################################################
import sys
sys.path.append("/opt/WPG/") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/WPG") # DESY MAXWELL PATH

sys.path.append("/opt/spb_model") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/spb_model") # DESY MAXWELL PATH
###############################################################################
###############################################################################
import numpy as np
import os
from wpg import srwlib
from wpg.wavefront import Wavefront
from wpg.wpg_uti_wf import integral_intensity, calculate_fwhm, look_at_q_space
from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity
from model.src.coherent import coherentSource
from model.beamline.structure import BeamlineModel, propParams
from wpg.srwlib import SRWLOptD as Drift
from wpg.srwlib import SRWLOptA as Aperture
from wpg.generators import build_gauss_wavefront
from wpg.beamline import Beamline

def loadWavefront(fname, indir = None):
    
    wfr = Wavefront()
    
    if indir is None:
        wfr.load_hdf5(fname)
    else:
        wfr.load_hdf5(indir + fname)
    
    return wfr
    

def constructPulse(nx = 512, ny = 512, nz = 512):
    
    wfr = Wavefront(build_gauss_wavefront(nx, ny, nz, 5.0, -400e-06, 400e-06, -400e-06, 400e-06, 1e-15, 5e-06, 5e-06, 19))
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
    #look_at_q_space(wfr)
    return wfr

def create_circular_mask(nx, ny, c = None, r = None):
    """
    create a circular mask numpy array
    
    :param nx: number of pixels of array (horizontal) [int]
    :param ny: number of pixels of array (vertical) [int]
    :param c: center of array in pixels [list]from wpg import srwlib
    :param r: radius of array (in pixels) [int]
    
    :returns mask: binary mask [np array]
    """ 
    if c is None: # use the middle of the image
        c = (int(nx/2), int(ny/2))
    if r is None: # use the smallest distance between the center and image walls
        r = min(c[0], c[1], nx-c[0], ny-c[1])

    X, Y = np.ogrid[:nx, :ny]
    dist_from_center = np.sqrt((X - c[0])**2 + (Y-c[1])**2)

    mask = dist_from_center <= r
    
    return mask


