# -*- coding: utf-8 -*-
import numpy as np
from felpy.model.beamlines.methods import get_beamline_object
from felpy.model.source.coherent import construct_SA1_pulse


if __name__ == '__main__':
    wfr = construct_SA1_pulse(1024, 1024, 5, 5.0, 0.25)
    spb = get_beamline_object(apertures = True, surface = "on")
    spb.propagate_sequential(wfr)