# -*- coding: utf-8 -*-

from wpg.generators import build_gauss_wavefront_xy as construct_gaussian
from felpy.model.core.wavefront import Wavefront
from wpg.wpg_uti_wf import plot_intensity_map 

from felpy.exp.NFS.cylinder_phase_mask import phase, cylinder_thickness
import numpy as np
from felpy.utils.opt_utils import ekev2wav
from felpy.model.tools import propagation_parameters
from wpg.optical_elements import Drift
from wpg.beamline import Beamline
from felpy.exp.NFS.speckle import define_speckle_mask

from felpy.analysis.statistics.correlation import norm


from PIL import Image


ekev = 25
wav = ekev2wav(25)
k = np.pi*2/wav
delta = 4e-07
d2waist = (wav)/(np.arctan(.27/8)*np.pi)

r = 75e-03




wfr = Wavefront(construct_gaussian(1024,1024, 25, -15e-02, 15e-02, -15e-02, 15e-02, 3e-2, 3e-02, d2waist))


### PLOT 1
plot_intensity_map(wfr)

bl = Beamline()

bl.append(Drift(1), propagation_parameters(1,1,1,1))
bl.propagate(wfr)


x = np.linspace(wfr.params.Mesh.xMin, wfr.params.Mesh.xMax, wfr.params.Mesh.nx)
y = np.ones([wfr.params.Mesh.nx, wfr.params.Mesh.ny])

thickness = norm(cylinder_thickness(x, r, offset = 0), lim = (0,255))

sdir = "../../data/samples/cylinder_thickness.png"

im = Image.fromarray(thickness*y)
im = im.convert("L")
im.save(sdir)

np.save(sdir, thickness)
rx, ry = wfr.get_spatial_resolution()

cylinder = define_speckle_mask(sdir, rx = rx, ry = ry,
                              sample_thickness = 2*r, sample_delta = 4e-07, plot = True, ekev = 25)


# =============================================================================
# phi_cyl = phase(x, r, k, delta, offset = -r)
# phi_cyl[np.isnan(phi_cyl)] = 0
# psi_cyl = 1+np.exp(1j*phi_cyl*y) #* np.exp(k*3e-10*hollow_cylinder_thickness(x, r))
# a = wfr.data.arrEhor[:,:,0,1]
# wfr.multiply_by_complex(psi_cyl)
# 
# =============================================================================
speckle = define_speckle_mask("../../data/samples/speckle_enlarged.tif", rx = 1.25e-06, ry = 1.25e-06,
                              sample_thickness = .85, sample_delta = 7e-08, plot = True, ekev = 25)


### PLOT 2
plot_intensity_map(wfr)
 
bl = Beamline()
#bl.append(cylinder, propagation_parameters(1/2,1,1/2,1))
bl.append(Drift(1), propagation_parameters(1,1,1,1, mode = 'quadratic'))
#bl.append(speckle, propagation_parameters(1/50,1,1/50,1))

bl.propagate(wfr)

### PLOT 3
plot_intensity_map(wfr,cmap = 'bone')

from matplotlib import pyplot as plt 


### PLOT 4
plot_intensity_map(wfr)

bl = Beamline()

bl.append(Drift(4.5), propagation_parameters(1,1,1,1, 'quadratic'))
bl.propagate(wfr)
plt.imshow(wfr.get_intensity().sum(-1), cmap = 'bone')
### PLOT 5
plot_intensity_map(wfr, cmap = 'bone')

from felpy.model.core.detector_test import Detector



### PLOT 6
plot_intensity_map(wfr)

d = Detector(6.5e-06, 6.5e-06, 2500,2500)
d.detect(wfr)

### PLOT 7
plot_intensity_map(wfr, cmap = 'bone')
