# -*- coding: utf-8 -*-
import numpy as np
from math import floor

from numpy.fft import fft2, fftshift, ifft2, ifftshift

def paganin_method(I, n_ratio, wav, z):
       
    Nx, Ny = I.shape
    dkx = 2 * np.pi / (Nx)
    dky = 2 * np.pi / (Ny)
    # 
    kx, ky = np.meshgrid((np.arange(0, Ny) - floor(Ny / 2) - 1) * dky, (np.arange(0, Nx) - floor(Nx / 2) - 1) * dkx)
    
    denom =  fft2(fftshift(I))
    A = (wav*z * n_ratio)/(4*np.pi)
    B = kx**2 + ky**2
    numer = 1/(1+A*B)
    numer[np.isnan(numer)] = 0
    
    
    return n_ratio/2  * np.log(ifftshift(ifft2(denom*numer)))

def test_paganin_method():
    
    I = load_tif("../../data/zyla_speckle/zyla_speckle_1.tif")
    n_ratio = 100
    wav = 0.5e-10
    z = 1.5 
    
    object_phase = paganin_method(I, n_ratio, wav, z)
    
    return I, object_phase

if __name__ == '__main__':
    
    from matplotlib import pyplot as plt
    from felpy.utils.np_utils import load_tif
    
    I, object_phase = test_paganin_method()
    plt.imshow(I.real, cmap = 'bone')
    plt.show()
    
    plt.imshow(object_phase.real, cmap = 'bone')
    plt.show()