# -*- coding: utf-8 -*-


import sys

from felpy.model.wavefront import Wavefront
from felpy.model.beamline import Beamline
from wpg.optical_elements import Drift
from felpy.model.tools import propagation_parameters
from felpy.utils.os_utils import mkdir_p


if __name__ == '__main__':
    print("working")
    in_directory = sys.argv[1]
    print("in: ", in_directory)
    out_directory = sys.argv[2]
    print("out: ", out_directory)
    wfr = Wavefront()
    wfr.load_hdf5(in_directory)
    print("wfr loaded")
    bl = Beamline()
    bl.append(Drift(3.644-2.2), propagation_parameters(1, 1, 1, 1, 'quadratic'))
    bl.propagate(wfr)
    print("wfr propagated")
    wfr.store_hdf5(out_directory)
    print("wfr stored")
    wfr.analysis(VERBOSE = True, DEBUG = True) 
    #print(wfr.custom_fields) ## good debug
        
    