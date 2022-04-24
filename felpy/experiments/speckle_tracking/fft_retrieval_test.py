# -*- coding: utf-8 -*-

import numpy as np
from numpy.fft import fft2, fftshift, ifft2, ifftshift

from math import floor as floor

# =============================================================================
# z2 = 1
# wav = 0.5e-10
# k = np.pi*2/wav
# 
# image = np.random.rand(500,500)
# dx = np.ones([500,500])
# dy = np.ones([500,500])*2
# Nx, Ny = image.shape
# dkx = 2 * np.pi / (Nx)
# dky = 2 * np.pi / (Ny)
# 
# kx, ky = np.meshgrid((np.arange(0, Ny) - floor(Ny / 2) - 1) * dky, (np.arange(0, Nx) - floor(Nx / 2) - 1) * dkx)
# ftfilt = 1 / (1j * kx - ky)
# ftfilt[np.isnan(ftfilt)] = 0
# C = (z2/k)
# 
# numerator = fft2(fftshift(dx + 1j*dy))
# delta_phi = C * ifftshift(ifft2((numerator)*(ftfilt)))
# =============================================================================
delta = 100
beta = 1 
a = delta/beta ### delta/beta
C = 0.5 * a
z2 = 1
wav = 0.5e-10
I0 = np.ones([500,500])
I1 = np.random.rand(500,500)


Nx, Ny = I0.shape
dkx = 2 * np.pi / (Nx)
dky = 2 * np.pi / (Ny)
# 
kx, ky = np.meshgrid((np.arange(0, Ny) - floor(Ny / 2) - 1) * dky, (np.arange(0, Nx) - floor(Nx / 2) - 1) * dkx)

denom =  fft2(fftshift(I1))
A = (wav*z2*delta)/(4*np.pi*beta)
B = kx**2 + ky**2
numer = 1/(1+A*B)
numer[np.isnan(numer)] = 0


ps = 1/2  * np.log(ifftshift(ifft2(denom*numer)))