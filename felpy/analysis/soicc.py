### second order intensity-intensity correlations

import numpy as np
from matplotlib import pyplot as plt

def calculate_coherence_length(g2, px = 1, VERBOSE = False):
    """
    calculate the coherence length via: 
    Kim, Y. et al.Journal of Synchrotron Radiation 29, no. 6 (November 1, 2022): 1465–79. 
    https://doi.org/10.1107/S1600577522008773.
    
    :param g2: two dimensional spatial coherence function g^{(2)}(r_1,r_2) of an ensemble of intensity distribution
    :param px: pixel size of the detector
    
    :returns L: coherence length in m
    """

    dx = g2.shape[0]//2
    g2de =np.diag(np.fliplr(g2))[dx:]
    delta = np.linspace(0,dx*px,dx)

    num = 2*np.sum([delta[itr]**2*(g2de[itr]-1) for itr in range(g2de.shape[0])])
    den = np.sum([(g2de[itr]-1) for itr in range(g2de.shape[0])])
    L = np.sqrt(num/den)
    
    if VERBOSE:
        print("Coherence Len", L *1e3, "mm")
    
    return L


def calculate_contrast(g2, VERBOSE = False):
    """
    calculate the contrast function due to polychromaticity.
    
    :param g2: two dimensional spatial coherence function g^{(2)}(r_1,r_2) of an ensemble of intensity distribution
    
    Vartanyants, I. A., and A. Singer.
    New Journal of Physics 12, no. 3 (March 31, 2010): 035004. 
    https://doi.org/10.1088/1367-2630/12/3/035004.

    """
    contrast = ((g2[g2.shape[0]//2, g2.shape[1]//2]-1))
    
    if VERBOSE:
        print("Contrast: ", contrast)
    return contrast

def calculate_doc(g2, I, VERBOSE = False):
    """
    calculate degree of coherence of array from g2 function and input intensity
    
    :param g2: two dimensional spatial coherence function g^{(2)}(r_1,r_2) of an ensemble of intensity distribution
    :param I: intensity profile which g2 was calculated for
    """
    num = 0
    for x1 in range(g2.shape[0]):
        for x2 in range(g2.shape[1]):
            num += (g2[x1,x2]-1)*I.mean(0)[x1]*I.mean(0)[x2]


    den = np.sum(I.mean(0))**2
    contrast = calculate_contrast(g2)
    doc = 1/contrast*(num/den)*100

    if VERBOSE:
        print("degree of coherence (%):", doc)
    
    return doc

def calculate_g2(i1,i2):
    """
    gets the 2d second-order intensity-intensity spatial correlation (soiic) of an electric field in time
    
    :param i1: intensity distribution 1 - [nt,nr]
    :param i2: intensity distribution 2 - [nt, nr]
 
    :returns g2: soicc - [x_1, x_2]
    """
 

    assert i1.shape == i2.shape, "Arrays should be of the same dimensions"
    
    idx = np.arange(0, i1.shape[1]) ### index over spatial domain (nr)

        
    g2 = np.zeros([idx.shape[0], idx.shape[0]]) ### init array
    mid = g2.shape[0]//2 ## mid point of time array

    for x1 in idx:

        r1 = i1[:, x1]
 
        for x2 in idx:


            r2 = i1[:, x2]

            num = np.mean(r1*r2)
            den = np.mean(r1)*np.mean(r2)
            corr = np.mean(num/den)
            g2[x1,x2] = corr
    
    return g2 


### plotters
def plot_g2(g2, axes = None):
    """
    plot the g2(x_1, x_2) and g2(\delta x) functions
    
    :param g2: two dimensional spatial coherence function g^{(2)}(r_1,r_2) of an ensemble of intensity distribution
    :param axes: [optional] provide a set of axes for plotting, list of len 2 mpl.axes objects
    """
    if axes is None:
        fig, axes = plt.subplots(1,2)
    
    else:
        assert type(axes) == list
        assert len(axes) == 2

    ax1, ax2 = axes
    
    ax1.imshow(np.fliplr(g2), vmin = 1)
 
     
    ax2.plot(np.diag(np.fliplr(g2)))
    
def gaussian_a(x,a,mu, sig):
    return 1+a*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))