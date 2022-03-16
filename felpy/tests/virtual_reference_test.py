# -*- coding: utf-8 -*-

from felpy.model.source.coherent import construct_SA1_wavefront
from wpg.srwl_uti_smp import srwl_opt_setup_transm_from_file as srwloptT
from felpy.utils.opt_utils import get_required_distance, get_magnification, fresnel_criterion, ekev2wav

import numpy as np
from matplotlib import pyplot as plt
from wpg.wpg_uti_wf import plot_intensity_map as plot

def odd_even(N):
    """
    check if a number is odd or even
    """
    
    if N % 2 == 0:
        ans = "even"
    elif N % 2 == 1:
        ans = 'odd'
    return ans

def shift_image(X, dx, dy):
   X = np.roll(X, dy, axis=0)
   X = np.roll(X, dx, axis=1)
   if dy>0:
       X[:dy, :] = 0
   elif dy<0:
       X[dy:, :] = 0
   if dx>0:
       X[:, :dx] = 0
   elif dx<0:
       X[:, dx:] = 0
   return X

def get_window(arr, c, l):
    """
    return the window of analysis given a feature pixel centered at C
    
    :param c: coordinates (tuple of window center) (cx,cy)
    :param l: side length of square window
    """
    cx = c[1]
    cy = c[0]
    
    pad = l//2
    f = np.pad(arr, l)
    fl = l//2
    
    if odd_even(l) == 'odd':
        feature = f[(cx+pad-l):(cx+pad),
                    (cy+pad):(cy+pad+l)]
 
    elif odd_even(l) == 'even':
        feature = f[cx+pad-l//2:cx+pad+l//2,
                    cy+pad-l//2:cy+pad+l//2]
    return feature



def window_stack(a, stepsize=1, width=3):
    n = a.shape[0]
    return np.hstack( [a[i:1+n+i-width:stepsize] for i in range(0,width)] )

wfr = construct_SA1_wavefront(512, 512, 5.0, q = 0.1, mx = 0, my = 0)
w = 15e-06
wav = wfr.get_wavelength()
px = wfr.pixelsize()
z = get_required_distance(w, wfr.get_spatial_resolution()[0], wav)
print(z)
f = fresnel_criterion(w, z*2, wav)
print(f)
z *= 2
 
from felpy.model.beamline import Beamline
from wpg.optical_elements import Drift
from felpy.model.tools import propagation_parameters

obj = srwloptT("/opt/FELpy/felpy/data/samples/lena.png", 0.33e-06, 0.33e-06, 20e-06, 3.3269971E-05, 60e-06)
bl = Beamline()
#bl.append(Drift(2), propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))

bl.append(obj, propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
bl.append(Drift(z), propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))

bl.propagate(wfr)
plot(wfr)
arr1 =  wfr.get_intensity() 

wfr = construct_SA1_wavefront(512, 512, 5, q = 0.1, mx = 0, my = 0, tiltX = 1e-06)
bl = Beamline()
plot(wfr)

#bl.append(Drift(2), propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))

bl.append(obj, propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
bl.append(Drift(z), propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))

bl.propagate(wfr)
arr2 =  wfr.get_intensity() 

arr1 = arr1[200:400,200:400,0]
arr2 = arr2[200:400,200:400,0]

#arr2 = shift_image(arr1, 2, 2)
plt.imshow(arr1[:,:])
# =============================================================================
# plt.show()
# plt.imshow(arr2[:,:])
# plt.show()
# 
# 
# c = np.zeros_like(arr1)
# from skimage.registration import phase_cross_correlation
# 
# for nx in range(arr1.shape[0]):
#     for ny in range(arr1.shape[1]):
#         #print(phase_cross_correlation(get_window(arr1, (nx,ny), 10), get_window(arr2, (nx,ny), 10)))
#         c[nx,ny] = phase_cross_correlation(get_window(arr1, (nx,ny), 10), get_window(arr2, (nx,ny), 10))[0][1]
#         #print(c[nx,ny])
#         plt.imshow(c)
# =============================================================================
