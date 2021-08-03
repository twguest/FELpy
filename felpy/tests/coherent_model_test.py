# -*- coding: utf-8 -*-
import numpy as np
from felpy.model.beamlines.exfel_spb.methods import get_beamline_object
from felpy.model.src.coherent import construct_SA1_wavefront, construct_SA1_pulse
from felpy.model.tools import propagation_parameters
from wpg.wpg_uti_wf import plot_intensity_map
if __name__ == '__main__':
    
    wfr = construct_SA1_wavefront(1024, 1024, 5.0, 0.25)
    
    spb = get_beamline_object(apertures = True, surface = "on", crop = ["d1","NVE"])
    spb.propagation_options[0]['optical_elements'][-1].L    += 2.2
    spb.propagation_options[0]['propagation_parameters'][-1] = propagation_parameters(1,1,1,1,mode = 'quadratic')
    spb.propagate(wfr)

    plot_intensity_map(wfr)
    
# =============================================================================
#     wfr = construct_SA1_pulse(2028, 2028,5, 10.0, 1.5)
#     
#     spb = get_beamline_object(apertures = True, surface = "on", crop = ["d1","NVE"])
#     spb.propagation_options[0]['optical_elements'][-1].L    += 2.2
#     spb.propagation_options[0]['propagation_parameters'][-1] = propagation_parameters(1,1,1,1,mode = 'quadratic')
#     spb.propagate_sequential(wfr)
# 
#     plot_intensity_map(wfr)
# =============================================================================
