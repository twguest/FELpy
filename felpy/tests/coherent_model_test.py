# -*- coding: utf-8 -*-
import numpy as np
from felpy.model.beamlines.exfel_spb.methods import get_beamline_object
from felpy.model.src.coherent import construct_SA1_wavefront


if __name__ == '__main__':
    wfr = construct_SA1_wavefront(1024, 1024, 5.0, 0.25)
    
    spb = get_beamline_object(apertures = True, surface = "on")
    spb.propagate_sequential(wfr)