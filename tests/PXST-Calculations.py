# -*- coding: utf-8 -*-
"""
@author: twguest

Rough scratchpad for PXST calculations
"""

import numpy as np

def magnification(src2sample, sample2det):
    
    M = (src2sample+sample2det)/src2sample
    
    return M

def phaseSensitivity(sigma_det, z, M, wav):
    """
    calculate effective phase sensitivity
    """
    ps = (np.pi*2)/wav * (sigma_det**2/(z*M))
    return ps

def effectiveResolution(interpolation, sigma_det, feature_size, z, z1, focalLength):
    """
    calculate effective resolution
    """
    
    sigma_eff = interpolation*(sigma_det + feature_size*((z*focalLength)/z1))
    return sigma_eff

def main():
    
    npixels = [2560, 2160]
    
    ### taking effective pixel resolution via Richard
    eff_x = 6.5e-06 
    eff_y = 6.5e-06 
    print("Effective Pixel Size (x/y): {:2e}, {:2e}".format(eff_x, eff_y))

    wav = 2.38e-10
    
    z1 = 0.01 # source to sample
    print("Source to Sample Distance: {}".format(z1))
    
    z = 3.5 # sample to det
    print("Sample to Detector Distance: {}".format(z))
    
    M = magnification(z1, z)
    print("Magnification: {}".format(M))
    
    ### COMPARATIVELY
    ## via simulation
    pix = 5e-10 ## approx effective pixel size of sample-plane wavefield
    print("Alternative Magnification: {}".format(eff_x/pix))
    
    
    fx = 3.2 #focal length of horizontal kb
    fy = 2.2 #focal length of vertical kb
    
    interpolation = 1
    print("Interpolation Factor: {}".format(interpolation))
    
    feature_size = 10e-09
    print("Feature Size: {:.2e} m".format(feature_size))
    

    psx = phaseSensitivity(eff_x, z, M, wav)
    psy = phaseSensitivity(eff_y, z, M, wav)
    
    print("Hor. Phase Sensitivity: {:.2e} rad".format(psx))
    print("Ver. Phase Sensitivity: {:.2e} rad".format(psy))
    
    refx = psx/M
    refy = psy/M
    
    print("Hor. Reference Resolution: {:.2e} m".format(refx))
    print("Ver. Reference Resolution: {:.2e} m".format(refy))
    
    angSens = [eff_x/(z*M), eff_y/(z*M)]
    print("Angular Resolution: {:.2e}, {:2e}".format(angSens[0], angSens[1]))
    

if __name__ == '__main__':
    main()