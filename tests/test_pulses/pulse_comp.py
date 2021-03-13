#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 11:48:50 2020

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

from felpy.model.core.wavefront import Wavefront
from wpg.wpg_uti_wf import integral_intensity, calculate_fwhm

from model.src.coherent import construct_SA1_wavefront
from model.beamline.structure import Instrument

def loadPulse(wdir):
    
    wfr = Wavefront()
    wfr.load_hdf5(wdir)
    
    wfr.initial_intensity = integral_intensity(wfr, bPlot = False)
    return wfr

def propagatePulse(wfr, outdir, mode = 'direct'):
    
    spb = Instrument()
    
    spb.setupHOMs(wfr.params.photonEnergy/1000, 2.2e-03)
    spb.setupKBs(wfr.params.photonEnergy/1000, 3.5e-03)
    
    spb.mirrorProfiles(toggle = "on", aperture = True, overwrite = True)
    
    spb.build_elements(focus = "micron")
<<<<<<< HEAD
    spb.build_beamline(focus = "micron")
=======
    spb.buildBeamline(focus = "micron")
>>>>>>> 108cfb9b6fc97d3841ee1db54862523eee5b184e
    
    #spb.scale(wfr, isc = 501) 
    
    bl = spb.get_beamline()
    
    if mode == 'direct':
        bl.propagate(wfr)
    elif mode == 'sequential':
        bl.propagate_sequential(wfr, outdir)
    
    wfr.final_intensity = integral_intensity(wfr, bPlot = False)
    
    return wfr

def constructCoherentEquiv(wfr, outdir, mode = 'direct'):
    
    cwfr = construct_SA1_wavefront(1024, 1024, wfr.params.photoEnergy/1000, 0.25)
    
    cwfr.initial_intensity = integral_intensity(cwfr)
    
    spb = Instrument()
    
    spb.setupHOMs(wfr.params.photonEnergy/1000, 2.2e-03)
    spb.setupKBs(wfr.params.photonEnergy/1000, 3.5e-03)
    
    spb.mirrorProfiles(toggle = "on", aperture = True, overwrite = True)
    
    spb.build_elements(focus = "micron")
<<<<<<< HEAD
    spb.build_beamline(focus = "micron")
=======
    spb.buildBeamline(focus = "micron")
>>>>>>> 108cfb9b6fc97d3841ee1db54862523eee5b184e
    
    #spb.scale(cwfr, isc = 150) 
    
    bl = spb.get_beamline()
    
    if mode == 'direct':
        bl.propagate(cwfr)
    elif mode == 'sequential':
        bl.propagate_sequential(cwfr, outdir)
    
    return cwfr

def comparePulses(wfr, cwfr, savedir):
    
    print("Intial Incoherent Integrated Intensity: {} W/mm^2".format(wfr.initial_intensity))
    print("Coherent Integrated Intensity: {} W/mm^2".format(cwfr.initial_intensity))
    
    print("Final Incoherent Integrated Intensity: {} W/mm^2".format(wfr.final_intensity))
    print("Final Coherent Integrated Intensity: {} W/mm^2".format(cwfr.final_intensity))
    
    print("Incoherent System Efficiency: {}".format((wfr.final_intensity)/wfr.initial_intensity))
    print("Coherent System Efficiency: {}".format((cwfr.final_intensity)/cwfr.initial_intensity))
    
    
    fwhm_x, fwhm_y = calculate_fwhm(wfr)['fwhm_x'], calculate_fwhm(wfr)['fwhm_y']
    cfwhm_x, cfwhm_y = calculate_fwhm(cwfr)['fwhm_x'], calculate_fwhm(cwfr)['fwhm_y']
    
    print("Incoherent Focal FWHM-x: {} um".format(fwhm_x*1e6))
    print("Incoherent Focal FWHM-y: {} um".format(fwhm_y*1e6))
    
    print("Coherent Focal FWHM-x: {} um".format(cfwhm_x*1e6))
    print("coherent Focal FWHM-y: {} um".format(cfwhm_y*1e6))
    
    ### store pulses
    wfr.store_hdf5(savedir + "pulse.hdf5")
    cwfr.store_hdf5(savedir + "gaussian.hdf5")
    
if __name__ == "__main__":
    
    wdir = "../../data/h5/gauss.h5"
    
    wfr = loadPulse(wdir)
    print(wfr.data.arrEhor.shape)
    #wfr = propagatePulse(wfr, outdir = "/tmp/", mode = 'sequential')
    
    #cwfr = constructCoherentEquiv(wfr, outdir = "/tmp/")
    
    #comparePulses(wfr, cwfr, "/data/")