# -*- coding: utf-8 -*-

import numpy as np
import seaborn as sns

from felpy.utils.np_utils import get_mesh
from felpy.utils.vis_utils import basic_plot, simple_line_plot, colorbar_plot, signal_plot, scatter_plot
from wpg.beamline import Beamline
from wpg.optical_elements import Aperture, Drift
from felpy.model.core.wavefront import Wavefront

from felpy.model.src.coherent import construct_SA1_wavefront

from matplotlib import pyplot as plt

from numpy import fft
from felpy.utils.opt_utils import ekev2wav

from scipy.constants import c





def wavefront_from_array(cfr,nx,ny,nz,dx,dy,dz,ekev, pulse_duration = 40e-15, sigma = 4):


    # Initialize empty wavefront.
    wfr = Wavefront()

    # Setup E-field.
    wfr.data.arrEhor = np.zeros(shape=(nx, ny, nz, 2))
    wfr.data.arrEver = np.zeros(shape=(nx, ny, nz, 2))

    wfr.params.wEFieldUnit = 'sqrt(W/mm^2)'
    wfr.params.photonEnergy = ekev * 1000
    wfr.params.wDomain = 'time'
    wfr.params.Mesh.nSlices = nz
    wfr.params.Mesh.nx = nx
    wfr.params.Mesh.ny = ny


    
    wfr.params.Mesh.sliceMin = -pulse_duration*sigma / 2.
    wfr.params.Mesh.sliceMax = pulse_duration*sigma / 2.

    range_x = dx*nx
    range_y = dy*ny

    wfr.params.Mesh.xMin = -range_x / 2.
    wfr.params.Mesh.xMax = range_x / 2.
    wfr.params.Mesh.yMin = -range_y / 2.
    wfr.params.Mesh.yMax = range_y / 2.


    wfr.data.arrEhor = complex_to_wpg(cfr)
    
    #wfr.set_electric_field_representation('f')
        
    return wfr

def wavefront_to_wavefield(spatial_wfr, temporal_profile):
    
    new_wfr = spatial_wfr.as_complex_array()[:,:,:]*temporal_profile


    dx, dy = spatial_wfr.get_spatial_resolution()
    dz = spatial_wfr.get_temporal_resolution()
    

    wfr = wavefront_from_array(new_wfr,
                             nx = spatial_wfr.params.Mesh.nx,
                             ny = spatial_wfr.params.Mesh.ny,
                             nz = len(temporal_profile),
                             dx = dx,
                             dy = dy,
                             dz = dz,
                             ekev = spatial_wfr.params.photonEnergy/1000)
    return wfr


def complex_to_wpg(arr): ### converter
    new_arr = np.zeros([arr.shape[0], arr.shape[1], arr.shape[2], 2])
    new_arr[:,:,:,0] = arr.real
    new_arr[:,:,:,1] = arr.imag
    return new_arr


if __name__ == '__main__':

    pulse_time = 100e-15
    sigma = 4 
    S = 4
    Seff = S/sigma
    
    n_samples, sampling_interval_t = temporal_sampling_requirements(pulse_time, VERBOSE = True, S = S)
    sampling_interval_w = 1/sampling_interval_t

    n_samples *= sigma 
    
    t = np.arange(-pulse_time*4, pulse_time*4, pulse_time/n_samples)

    temporal_profile = generate_temporal_SASE_pulse(pulse_time = pulse_time,
                                                    n_samples = n_samples,
                                                    sigma = sigma,
                                                    VERBOSE = True)
    print("Number of Samples: ",n_samples)
    print(temporal_profile.dtype)
    spatial_profile = construct_SA1_wavefront(128, 128, 5.0, 0.25)

    wfr = wavefront_to_wavefield(spatial_profile, temporal_profile)
    
    from wpg.wpg_uti_wf import plot_intensity_map
    plot_intensity_map(wfr)
    wfr.set_electric_field_representation('f')
    from felpy.model.beamlines.exfel_spb.methods import setup_spb
    spb = setup_spb(parameter_file = "/opt/FELpy/felpy/data/params/spb-sfx_nkb_FAST.json", theta_KB = 5e-03, theta_HOM = 3.5e-03) #bl = spb.bl
    bl = spb.bl 
    
    ## EXEC
    bl.propagate_sequential(wfr)
