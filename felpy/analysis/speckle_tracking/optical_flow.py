import glob
import os
import sys
from numpy.fft import fftshift as fftshift
from numpy.fft import ifftshift as ifftshift
from numpy.fft import fft2 as fft2
from numpy.fft import ifft2 as ifft2
from numpy.fft import fftfreq as fftfreq

from scipy.ndimage.filters import gaussian_filter
from math import pi as pi
from math import floor as floor

from felpy.analysis.speckle_tracking import frankoChellappa as fc



import numpy as np
 
def optical_flow(intensityImage,derivative,alpha=0,sig_scale=0):

    Nx, Ny = derivative.shape
    # fourier transfomm of the derivative and shift low frequencies to the centre
    ftdI = fftshift(fft2(derivative))
    # calculate frequencies
    dqx = 2 * pi / (Nx)
    dqy = 2 * pi / (Ny)
    Qx, Qy = np.meshgrid((np.arange(0, Ny) - floor(Ny / 2) - 1) * dqy, (np.arange(0, Nx) - floor(Nx / 2) - 1) * dqx)


    #building filters


    sigmaX = dqx / 1. * np.power(sig_scale,2)
    sigmaY = dqy / 1. * np.power(sig_scale,2)
    #sigmaX=sig_scale
    #sigmaY = sig_scale

    g = np.exp(-(((Qx)**2) / 2. / sigmaX + ((Qy)**2) / 2. / sigmaY))
    #g = np.exp(-(((np.power(Qx, 2)) / 2) / sigmaX + ((np.power(Qy, 2)) / 2) / sigmaY))
    beta = 1 - g;

    # fourier filters
    ftfiltX = (1j * Qx / ((Qx**2 + Qy**2 + alpha))*beta)
    ftfiltX[np.isnan(ftfiltX)] = 0

    ftfiltY = (1j* Qy/ ((Qx**2 + Qy**2 + alpha))*beta)
    ftfiltY[np.isnan(ftfiltY)] = 0

    # output calculation
    dImX = 1. / intensityImage * ifft2(ifftshift(ftfiltX * ftdI))
    dImY = 1. / intensityImage * ifft2(ifftshift(ftfiltY * ftdI))

    return dImX.real,dImY.real




def kottler(dX,dY):
    print('kottler')
    i = complex(0, 1)
    Nx, Ny = dX.shape
    dqx = 2 * pi / (Nx)
    dqy = 2 * pi / (Ny)
    Qx, Qy = np.meshgrid((np.arange(0, Ny) - floor(Ny / 2) - 1) * dqy, (np.arange(0, Nx) - floor(Nx / 2) - 1) * dqx)

    polarAngle = np.arctan2(Qx, Qy)
    ftphi = fftshift(fft2(dX + i * dY))*np.exp(i*polarAngle)
    ftphi[np.isnan(ftphi)] = 0
    phi3 = ifft2(ifftshift(ftphi))
    return phi3.real



def LarkinAnissonSheppard(dx,dy,alpha =0 ,sigma=0):
    Nx, Ny = dx.shape
    i = complex(0, 1)
    G= dx + i*dy
    # fourier transfomm of the G function
    fourrierOfG = fftshift(fft2(G))


    dqx = 2 * pi / (Nx)
    dqy = 2 * pi / (Ny)
    Qx, Qy = np.meshgrid((np.arange(0, Ny) - floor(Ny / 2) - 1) * dqy, (np.arange(0, Nx) - floor(Nx / 2) - 1) * dqx)

    ftfilt = 1 / (i * Qx - Qy)
    ftfilt[np.isnan(ftfilt)] = 0
    phi=ifft2(ifftshift(ftfilt*fourrierOfG))
    phi=phi.real
    return phi




def process_optical_flow(Is,Ir, sigma =1, alpha = .1):

    dI = (Is - Ir * (np.mean(gaussian_filter(Is,sigma=sigma)) / np.mean(gaussian_filter(Ir,sigma=sigma))))
    
    alpha=np.finfo(np.float32).eps
    dx, dy = optical_flow(Is, dI, alpha=alpha, sig_scale=sigma)
    
    return dx, dy

def process_all(Is,Ir, sigma =1, alpha = .1):

    dI = (Is - Ir * (np.mean(gaussian_filter(Is,sigma=sigma)) / np.mean(gaussian_filter(Ir,sigma=sigma))))
    
    alpha=np.finfo(np.float32).eps
    dx, dy = optical_flow(Is, dI, alpha=alpha, sig_scale=sigma)
    
    phi = fc.frankotchellappa(dx, dy, False)
    phi3 = kottler(dx, dy)
    phi2 = LarkinAnissonSheppard(dx, dy)


    return {'dx': dx, 'dy': dy, 'phi': phi, 'phi2': phi2,'phi3': phi3}



if __name__ == "__main__":
    pass