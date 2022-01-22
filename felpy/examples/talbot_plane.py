# -*- coding: utf-8 -*-

from felpy.model.source.coherent import construct_SA1_wavefront
from felpy.model.core.beamline import Beamline

from wpg.srwl_uti_smp import srwl_opt_setup_transm_from_file as srwloptT

from felpy.model.tools import create_circular_mask
from matplotlib import pyplot as plt
import numpy as np
from PIL import Image

from wpg.optical_elements import Drift
from felpy.model.tools import propagation_parameters

from wpg.wpg_uti_wf import plot_intensity_map 
def create_mask(nx, ny):
    """
    note nx and ny are the dimensions of a single mask element
    """
    arr = np.ones([nx, ny])*256
    arr[arr.shape[0]//2-5:arr.shape[0]//2+5, arr.shape[1]//2-5:arr.shape[1]//2+5] = 0
    arr = np.tile(arr, [5,5])
    return arr
    
nx = 25
ny = 25

plt.imshow(create_mask(nx, ny))
wfr = construct_SA1_wavefront(512, 512, 5, .8, mx = 0, xoff = 0, tiltX = 0)


period = nx * wfr.get_spatial_resolution()[0]*(512/(nx))
print(period)
wav = wfr.get_wavelength()

d = (2*period**2)/wav

arr = create_mask(nx, ny)
np.save("/opt/arr", arr)
im = Image.fromarray(arr)
im = im.convert('RGB')
im.save("/opt/arr.png")

obj = srwloptT("/opt/arr.png", 
               wfr.get_spatial_resolution()[0]*(512/75),
               wfr.get_spatial_resolution()[1]*(512/75),
               1e-06,
               1e-8,
               atten_len = 1e-08)

bl = Beamline()
bl.append(Drift(50),propagation_parameters(1,1,1,1, mode = 'quadratic') )
bl.append(obj, propagation_parameters(1, 1,1,1))
bl.append(Drift(d),propagation_parameters(1/3,4,1/3,4, mode = 'fraunhofer') )
bl.propagate(wfr)
plot_intensity_map(wfr)