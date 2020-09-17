#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 05:52:37 2020

@author: guestt
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""2
Created on Tue Jun 30 15:09:54 2020

Python Operations for NanoKB test w/ FAST PULSE
@author: twguest
"""

import sys

sys.path.append("../../")
import os
 
from tqdm import tqdm

from model.beamline.structure import BeamlineModel, propParams
from model.src.coherent import coherentSource
 
from wpg.beamline import Beamline
from wpg.srwlib import SRWLOptD as Drift

from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity

from model.h5_tools import obj, object2h5, h52object  
from wpg.multisliceOptE import greyscaleToSlices

import numpy as np 
from utils.job_utils import JobScheduler
from model.materials.phaseMask import phaseMask
from scipy.ndimage import gaussian_filter
from utils.os_utils import mkdir_p

from wpg.srwl_uti_smp import srwl_opt_setup_transm_from_file as Sample

def phaseScreen(wfr):
    
    
    wav = 2e-10
    ps = np.random.rand(250, 250)*wav
    ps *= gaussian_filter(ps, sigma = 75)
    
    opd = phaseMask(ps, [20e-06, 20e-06], wav)
    
    bl = Beamline()
    bl.append(opd, propParams(1,1,1,1, mode = 'normal'))
    bl.propagate(wfr)
    
    return wfr
    

def propThruMask(wfr, _x = 0, _y = 0):
    
    slices = greyscaleToSlices("../../data/samples/AAO.png",
                               propagationParameters = propParams(1,1,1,1, mode = 'normal'),
                               slices = 1,
                               resolution = [10e-09, 10e-09,  None],
                               thickness = 60e-06,
                               delta = 3.35e-03, #### via Henke
                               attenLength = 20e-06,
                               rotate_angle = 0,
                               shift_x = _x,
                               shift_y = _y)
    
    slices.propagate(wfr)
    
def propThruMaskLite(wfr, _x = 0, _y = 0):
    """
    propagate through the speckle generator
    """
    s = Sample("../../data/samples/AAO.png", 
               rx = 1e-09, ry = 1e-09,
               thickness = 60e-06,
               delta = 3.35e-03,
               atten_len = 20e-06,
               xc = 0, yc = 0,
               shift_x = _y, shift_y = _x,
               tile = [5,5])


    bl = Beamline()
    bl.append(s, propParams(1,1,1,1,mode = 'normal'))
    bl.propagate(wfr)
    
    return wfr
    

def propagateToMask(wfr, foc2sample):
    
    """
    propagate from the focal plane to the mask
    
    :param wfr: wpg wavefront structure
    :param foc2sample: focus to sample distance [m] (float)
    """
    
    focus = "nano"
    
    spb = BeamlineModel()
    
    spb.setupHOMs(wfr.params.photonEnergy/1000, 2.2e-03)
    spb.setupKBs(wfr.params.photonEnergy/1000, 3.5e-03)
    
    spb.mirrorProfiles(toggle = "on", aperture = True, overwrite = False)
    
    spb.buildElements(focus)
    spb.buildBeamline(focus)
    print("Setup Beamline")
    

    bl = spb.get_beamline()
    bl.propagation_options[0]['optical_elements'][-1].L += foc2sample
    bl.propagation_options[0]['propagation_parameters'][-1] = propParams(1, 1, 1, 1, mode = 'quadratic')

    bl.propagate(wfr)
    
    return wfr

 




def propagate2detector(wfr, sample2detector):
    
    
  
    bl = Beamline()
    bl.append(Drift(sample2detector), propParams(2, 1, 1, 1, mode = 'quadratic'))
    bl.propagate(wfr)
    
    return wfr

def resizeImage(wfr, nx = None, ny = None, fov_x = None, fov_y = None):
    
    bl = Beamline()
    
    if nx == None:
        nx = 2560
    if ny == None:
        ny = 2160
    if fov_x == None:
        fov_x = nx*6.5e-06
    if fov_y == None:
        fov_y = ny*6.5e-06
    
    [xMin, xMax, yMin, yMax] = wfr.get_limits()
    bl.append(Drift(0), propParams(fov_x/(xMax-xMin), nx/wfr.params.Mesh.nx, fov_y/(yMax-yMin), ny/wfr.params.Mesh.ny, mode = 'normal'))
    bl.propagate(wfr)
    
    return wfr

def get_random_string(length):
    """ TEMPORARY LOCATION
    writes a random string with length characters
    """
    
    import string
    import random
    
    letters = string.ascii_letters
    result_str = ''.join(random.choice(letters) for i in range(length))
    return result_str


def launch():
    """
    This part launches the jobs that run in main 
    """
    
    cwd = os.getcwd()
    script = "PXSE_experiment.py"
    
    js = JobScheduler(cwd + "/" + script, logDir = "../../logs/",
                      jobName = "coherentPXST", partition = 'exfel', nodes = 1, jobType = 'spawn', nSpawn = 25)
    
        
    js.run(test = True)
    
    

def kirkwoodTest():
        
    px = 1e-09
    
    for i in tqdm(range(-5, 5, 1)):
        
        print("Translating Mask at 5 Micron Steps")

        pos_x = i*px
        pos_y = 0
                        
        mkdir_p("/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/PXST/")
        outdir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/PXST/kirkwood/"
        
        mkdir_p(outdir)
        
        wfr = coherentSource(2560, 2160, 4.96, 0.25)
        print("Wavefront Loaded")
        wfr = propagateToMask(wfr, 1000e-06)
        wfr.save_tif(outdir + "atMask_{}.tif".format(i))
        
        wfr = phaseScreen(wfr)
        
        propThruMask(wfr, _x = pos_x, _y = pos_y)
        wfr.save_tif(outdir + "postMask_{}.tif".format(i))
        
        sample2detector = 4.0
        wfr = propagate2detector(wfr, sample2detector)
        
    
        wfr = resizeImage(wfr)
        wfr.save_tif(outdir + "atDet_{}.tif".format(i))        


def testPXST(xpos = 0, ypos = 0, nx = 2560, ny = 2160,
             fov_x = None, fov_y = None,
             sam2foc = 1000e-06):

    print("Testing the coherent case of PXST")
    
    wfr = coherentSource(2560, 2160, 4.96, 0.25)
    print("Wavefront Loaded")
    
    wfr = propagateToMask(wfr, sam2foc)
    wfr = phaseScreen(wfr)
    
    propThruMaskLite(wfr, _x = xpos, _y = ypos)
    
    sample2detector = 4.0 - sam2foc
    wfr = propagate2detector(wfr, sample2detector)
    
    wfr = resizeImage(wfr,  nx = nx, ny = ny, fov_x = fov_x, fov_y = fov_y)
    plotIntensity(wfr)    
        
    
    

if __name__ == '__main__':
    
    import imageio
    
    seed = sys.argv[1] ##random seed
    np.random.seed = (seed)
    
    scanRange = 5e-06
    px = 10e-09
    
    delta = scanRange/px
    
    for i in range(15):
        
    
        pos_x = np.random.randint(-int(delta//2), int(delta//2))
        pos_y = np.random.randint(-int(delta//2), int(delta//2))
        
        data = obj()
        data.translation = [pos_x*px, pos_y*px]
        
        focus = "nano"
        mkdir_p("/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/PXST/")
        outdir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/PXST/coherentOut/"
        mkdir_p(outdir)
        
        wfr = coherentSource(2560, 2160, 4.96, 0.25)
        print("Wavefront Loaded")
        wfr = propagateToMask(wfr, 500e-06)
        
        wfr = phaseScreen(wfr)
        data.entryWavefront = wfr.toComplex()[0,:,:,0]
        data.x_pixel_size_i = wfr.pixelsize()[0]
        data.y_pixel_size_i = wfr.pixelsize()[1]
        
        
        propThruMask(wfr, _x = pos_x, _y = pos_y)
        
        sample2detector = 3.5
        wfr = propagate2detector(wfr, sample2detector)
        

        wfr = resizeImage(wfr)
        data.ComplexSolution = wfr.toComplex()[0,:,:,0]
        data.detectorIntensity = wfr.get_intensity().sum(axis = -1)
        
        lims = wfr.get_limits()
        data.x_pixel_size = wfr.pixelsize()[0]
        data.y_pixel_size = wfr.pixelsize()[1]
        data.distance = sample2detector
        data.mask = np.array(imageio.imread("../../data/samples/AAO.png"))
        data.wavelength = wfr.params.wavelength
        

        ID = get_random_string(10)
        object2h5(outdir + "SpeckleTest_{}{}.h5".format(seed, ID), data)
        ### DEBUG testpyobj = h52object(outdir + str(seed) + ".h5") twg 07/2020
        
        
    
