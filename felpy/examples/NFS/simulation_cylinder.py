# -*- coding: utf-8 -*-
 
from felpy.model.materials.phase_mask import phase_mask
import numpy as np
 
from felpy.utils.opt_utils import ekev2wav
from .examples.NFS.cylinder_phase_mask import phase 
from felpy.model.src.coherent import construct_gaussian

from wpg.wpg_uti_wf import plot_intensity_map as plot_ii
from wpg.optical_elements import Drift
from felpy.model.beamline import Beamline 
from felpy.model.tools import propagation_parameters
from .examples.NFS.speckle import define_speckle_mask

from felpy.model.detector_test import Detector

r1 = 10e-02
r2 = 5e-02

delta = 1e-04
ekev = 25
wav = ekev2wav(ekev)

k = (np.pi*2)/wav

x = np.linspace(-r1/2, r1/2, 100)

y = np.ones([100,100])

phase_shift = phase(x, r1, r2, k, delta)
mask = y*phase_shift



z1 = 2
z2 = .01

bl = Beamline()

from felpy.model.src.coherent import construct_SA1_wavefront

# =============================================================================
# src = construct_gaussian(1024, 1024, ekev, extent = [-25e-04, 25e-04, -25e-04, 25e-04],
#                      sigX = .5e-04, sigY = .5e-04, divergence = np.arctan(13.5e-03/8))
# 
# =============================================================================
q = 0.50
src = construct_SA1_wavefront(1024, 1024, 25.0, q)

plot_ii(src)

speckle = define_speckle_mask("../../data/samples/checkerboard-pattern.jpg", rx = 1e-8, ry = 1e-8,
                              sample_thickness = 150e-06, sample_delta = 4e-04, plot = True, ekev = 25)

bl.append(Drift(10), propagation_parameters(1,1,1,1, mode = 'quadratic'))

bl.append(speckle, propagation_parameters(1,1,1,1, mode = 'fresnel'))

#plot_ii(src)
bl.append(Drift(z2), propagation_parameters(1/4,4,1/4,4, mode = 'quadratic'))
bl.propagate(src)
src.view()
# =============================================================================
# plot_ii(src)
# 
# d = Detector(5e-07, 5e-07, 1024, 1024)
# d.detect(src)
# =============================================================================
plot_ii(src)
src.save_tif("/opt/test_{}nC".format(q))