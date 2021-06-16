# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "1.0.1"
__maintainer__ = "Trey Guest"
__email__ = "twguest@students.latrobe.edu.au"
__status__ = "Developement"

A script to plot analytical expectations of a phase contrast imaging experiment
using pyrex cylinder to introduce known phase gradients in a 25 keV x-ray beam
to be measured using NFS tracking techniques.
"""

import numpy as np
from matplotlib import pyplot as plt
from felpy.exp.NFS.cylinder_phase_mask import phase_gradient, phase
from felpy.utils.opt_utils import ekev2wav
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable


def fresnel_magnification(z1,z2):
    return (z1+z2)/z1

def angular_sensitivity(dx, dz):
    return np.arctan(dx/dz)





#ax1.set_xlim([0e-03, 15e-03])


def plot_fullfield_image(r):

    theta = phase_gradient(x, r, k, delta, offset = 0)
    theta[np.isnan(theta)] = 0 
    y = np.ones([1,N])*theta
        
    fig, ax1 = plt.subplots()
    plt.style.use(['science','ieee'])

    im1 = ax1.imshow(y*1e6,
                    extent = [x.min()*M*1e3, x.max()*M*1e3,
                              x.min()*M*1e3,x.max()*M*1e3],
                    aspect = 'auto', cmap = 'bone')
    ax1.set_yticks([])
    #ax1.set_xlim([-5, 15])
    
    ax1.set_xlabel("x (mm)")

    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes('right', size='7.5%', pad=0.05)
    
    cbar1 = fig.colorbar(im1, cax=cax1, orientation='vertical')
    cbar1.set_label("Phase Gradient ($\mu$rad)")

    ax1.annotate("$r$: {:2} mm".format(r*1e3), horizontalalignment = 'left',
             verticalalignment = 'bottom',
             xy = (-200,+200),
             c = 'black', fontsize = 10)
    
    
def plot_detector_image(r):
    theta = phase_gradient(x, r, k, delta, offset = r)
    theta[np.isnan(theta)] = 0 
    y = np.ones([1,N])*theta
    
    fig, ax1 = plt.subplots()
    
    im1 = ax1.imshow(y*1e6,
                    extent = [x.min()*M*1e3, x.max()*M*1e3,
                              x.min()*M*1e3,x.max()*M*1e3],
                    aspect = 'auto', cmap = 'bone', vmax = 0)
    
    ax1.set_yticks([])
    ax1.set_xlim([-15, 15])
    ax1.set_xlabel("x (mm)")

    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes('right', size='7.5%', pad=0.05)
    
    cbar1 = fig.colorbar(im1, cax=cax1, orientation='vertical')
    cbar1.set_label("Phase Gradient (mrad)")
 
    
    ax1.annotate("$r$: {:2} mm".format(r*1e3), horizontalalignment = 'left',
             verticalalignment = 'bottom',
             xy = (-15,200),
             c = 'black', fontsize = 10)
    
R = [1e-03, 5e-03, 10e-03, 25e-03, 50e-03, 75e-03]


### get fullfield images


def delta_phase_gradient(r1, r2):
    
    theta = abs(phase_gradient(x, r1, k, delta, offset = r1)-phase_gradient(x, r2, k, delta, offset = r2))
    theta[np.isnan(theta)] = 0 
    y = np.ones([1,N])*theta
    
    fig, ax1 = plt.subplots()
    
    im1 = ax1.imshow(y*1e6,
                    extent = [x.min()*M*1e3, x.max()*M*1e3,
                              x.min()*M*1e3,x.max()*M*1e3],
                    aspect = 'auto', cmap = 'jet')
    
    ax1.set_yticks([])
    ax1.set_xlim([-15, 15])
    ax1.set_xlabel("x (mm)")

    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes('right', size='7.5%', pad=0.05)
    
    cbar1 = fig.colorbar(im1, cax=cax1, orientation='vertical')
    cbar1.set_label("$\Delta\partial\phi$ ($\mu$rad)")
 
    
    ax1.annotate("$r_1$: {:2} mm".format(r1*1e3), horizontalalignment = 'left',
             verticalalignment = 'bottom',
             xy = (-15,200),
             c = 'black', fontsize = 10)
    ax1.annotate("$r_2$: {:2} mm".format(r2*1e3), horizontalalignment = 'left',
         verticalalignment = 'bottom',
         xy = (-15,160),
         c = 'w', fontsize = 10)

    



def plot_phase(r):
        
    phases = phase(x, r, k, delta, offset = r)
    phases[np.isnan(phases)] = 0 
    wrapped_phases = (np.pi * phases % (2 * np.pi) - np.pi)
    y = np.ones([1,N])*wrapped_phases
    fig, ax1 = plt.subplots()
    
    im1 = ax1.imshow(y,
                    extent = [x.min()*M*1e3, x.max()*M*1e3,
                              x.min()*M*1e3,x.max()*M*1e3],
                    aspect = 'auto', cmap = 'hsv', vmin = -np.pi, vmax = np.pi)
    
    
    ax1.set_yticks([])
    ax1.set_xlim([-15, 15])
    ax1.set_xlabel("x (mm)")

    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes('right', size='7.5%', pad=0.05)
    
    cbar1 = fig.colorbar(im1, cax=cax1, orientation='vertical')
    cbar1.set_label("$\Delta\partial\phi$ ($\mu$rad)")
    
    ax1.annotate("$r_1$: {:2} mm".format(r*1e3), horizontalalignment = 'left',
             verticalalignment = 'bottom',
             xy = (-15,200),
             c = 'black', fontsize = 10)
    
    return wrapped_phases



def plot_1d_phase_gradient(R):
    
    fig, ax1 = plt.subplots()
    
    for r in R:
        theta_r = phase_gradient(x, r, k, delta, offset = r)
        theta_r[np.isnan(theta_r)] = 0 

        ax1.plot(x*M*1e3, theta_r*1e6, label = "{} mm".format(r*1e3))
        
    ax1.set_xlim([-2.5,15])
    ax1.set_ylim([-30, 0])
    ax1.set_xlabel("x (mm)")
    ax1.set_ylabel("Phase Gradient ($\mu rad$")
    plt.legend()

def plot_1d_phase_difference(dr):
    
    fig, ax1 = plt.subplots()
     
    for r in dr:
        theta_r = phase_gradient(x, r, k, delta, offset = r)

        ax1.plot(x*M*1e3, theta_r*1e6, label = "{} mm".format(r*1e3))
        
    ax1.set_xlim([-2.5,15])
    ax1.set_ylim([-30, 0])
    ax1.set_xlabel("x (mm)")
    ax1.set_ylabel("Phase Gradient ($\mu rad$")
    plt.legend()




if __name__ == '__main__':
    delta = 4e-07 ### perspex
    ekev = 25
    wav = ekev2wav(ekev)
    
    k = (np.pi*2)/wav
    N = 2000
    
    x = np.linspace(-5e-02,5e-02,N)
    
    z1 = 2
    z2 = 6
    
    print("Angular Sensitivity: {}".format(angular_sensitivity(6.5e-06, z2)))
    
    M = fresnel_magnification(z1, z2)
    
    
    R = [15e-03, 25e-03, 50e-03, 75e-03]
    
    plot_fullfield_image(r = 15e-03)
    plot_detector_image(r = 15e-03)
    
    
    plot_fullfield_image(r = 75e-03)
    plot_detector_image(r = 75e-03)
    
    delta_phase_gradient(r1 = 15e-03, r2 = 75e-03)
    plot_1d_phase_gradient(R)