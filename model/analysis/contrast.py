#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 15:03:43 2020

@author: twguest


determine the contrast of an image by isolating spatial frequencies comparable
to the expected particle size
"""

import numpy as np
from matplotlib import pyplot as plt
import imageio 

import psychopy.visual
import psychopy.event
from copy import copy
from psychopy import filters
from matplotlib import pyplot as plt

def michelsonContrast(image):
    """
    determine the michelson contrast of an image
    
    :param image: image of interest (filepath or np-array)
    
    :returns contrast: contrast of image (float)
    """
    
    contrast = (np.max(image)-np.min(image))/(np.max(image)+np.min(image))

    return contrast

def bandpass(image, f1, f2, order = 10, bPlot = False):
    """
    bandpass filters to cut spatial frequencies from an image
    
    :param image: image to be filtered (filepath or np-array)
    :param f1: cut-in frequency [1/pix]
    :param f2: cut-off frequency [1/pix]
    :param order: order of the filter polynomial  
    
    :returns img_new: filtered image in r-space
    """
    
    if type(image) == str:
        image = imageio.imread(image)
    else:
        pass
    
    if bPlot:
        plt.imshow(image, cmap = 'bone')
        plt.title("Unfiltered Real-Space Image")
        plt.show()

    
    ### convert to frequency domain
    img_freq = np.fft.fft2(image)
    
    ### calculate amplitude spectrum
    img_amp = np.fft.fftshift(np.abs(img_freq))
    
    if bPlot:
        # for display, take the logarithm
        img_amp_disp = np.log(img_amp + 0.0001)
        
        plt.imshow(img_amp_disp, cmap = 'bone')
        plt.title("Unfiltered Frequency-Space Image")
        plt.show()
        
        
    ### construct filter
    f = filters.butter2d_bp(img_freq.shape, f1, f2, order)
    img_filt = np.fft.fftshift(img_freq) * f
    
        
    ### construct new image
    img_new = abs( np.real(np.fft.ifft2(np.fft.ifftshift(img_filt))) )
    
    
    
    if bPlot:
        img_filt = np.fft.fft2(img_new)
        img_filt = np.fft.fftshift(np.abs(img_filt))
        img_filt = np.log(img_filt + 0.0001)
        
        plt.imshow(img_filt.real, cmap = 'bone')
        plt.title("Filtered Frequency-Space Image")
        plt.show()
        
    if bPlot:
        plt.imshow(img_new, cmap = 'bone')
        plt.title("Filtered Real-Space Image")
        plt.show()
        
    return img_new    


def speckleContrast(image, speckleSize, ftol = 0.5, VERBOSE = True):
    """
    
    For determining the ideal plane of for a speckle mask
    
    :param image: detector image (filepath or np-array)
    :param speckleSize: expected speckle size (in pixels)
    :param ftol: filter tolerance (frac of spatial freq of speckles)
    :param VERBOSE: print to std.out (BOOL)
    
    :returns contrast: speckle contrast [0-1]
    """
    
    f1 = 1/speckleSize - ftol/speckleSize
    f2 = 1/speckleSize + ftol/speckleSize
    
    if VERBOSE: 
        bPlot = True
    else: 
        bPlot = False
        
    fim = bandpass(image, f1, f2, bPlot = bPlot)
    
    contrast = michelsonContrast(fim)
    return contrast


if __name__ == '__main__':
    
    
    ### TEST SPECKLE
    ### add increasing amounts of noise, should results in reduced contrast
    
    image = imageio.imread("../../data/samples/flower.png")  ### load an image    
    speckleSize = 10
    
    speckleContrast(image, speckleSize, VERBOSE = True)
