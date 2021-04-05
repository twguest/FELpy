# -*- coding: utf-8 -*-
import numpy as np
from felpy.model.beamlines.exfel_spb.methods import get_beamline_object
from felpy.model.src.coherent import construct_SA1_wavefront
from wpg.wpg_uti_wf import plot_intensity_map
from felpy.model.tools import propagation_parameters

if __name__ == '__main__':
    wfr = construct_SA1_wavefront(1024, 1024, 9.0, 0.25)
    spb = get_beamline_object(apertures = True, surface = "on")
    spb.propagation_options[0]['optical_elements'][-1].L += 2.2

    spb.propagation_options[0]['propagation_parameters'][-1] = propagation_parameters(2,1,2,1,mode = 'quadratic')
      
    spb.propagate_sequential(wfr)
 
    plot_intensity_map(wfr)