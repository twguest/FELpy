#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 14:17:32 2020

@author: twguest
"""

#################################################################
import sys
sys.path.append("/opt/WPG/") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/WPG") # DESY MAXWELL PATH

sys.path.append("/opt/spb_model") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/spb_model") # DESY MAXWELL PATH
###############################################################################
###############################################################################
from poppy import zernike
import numpy as np
from model.tools import constructPulse, create_circular_mask, mkdir_p
from felpy.model.core.wavefront import Wavefront

        
        
def fitZernicke(wfr, mode = 'integrated', nterms = 1250):
    
    nx, ny, nz = wfr.params.Mesh.nx, wfr.params.Mesh.ny, wfr.params.Mesh.nSlices
    
    aperture = create_circular_mask(nx, ny, r = nx//2-2)
    
    if mode == 'integrated':

        ph = wfr.data.arrEhor[:,:,:,1].sum(axis = 2)
        #ii = wfr.get_intensity()[:,:,0]

        zc = zernike.opd_expand(ph, aperture = aperture, nterms = nterms)
        #rc = zernike.opd_from_zernikes(zc)

        
    return zc

def test():
    wfr = constructPulse(nz = 10)
    zc = fitZernicke(wfr)
    return zc


if __name__ == '__main__':
    
    #ph, zc = test()
    
    mode = 'integrated'
    where = 'in' ## source or sample
    fname = sys.argv[1]
    
    
    outdir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/data/zernike/"
    savdir = outdir + mode + "/"
    
    mkdir_p(outdir)
    mkdir_p(savdir)
    
    wfr = Wavefront()
    wfr.load_hdf5("/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/{}/{}".format(where, fname))
    
    zc = fitZernicke(wfr, mode = mode, nterms = 1250)
    
    np.save(savdir + fname, zc)
