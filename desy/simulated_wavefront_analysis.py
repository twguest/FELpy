# -*- coding: utf-8 -*-


import sys

from felpy.model.core.wavefront import Wavefront

if __name__ == '__main__':
    print("working")
    in_directory = sys.argv[1]
    out_directory = sys.argv[2]
    
    wfr = Wavefront()
    print("wfr loaded")
    wfr.load_hdf5(in_directory)
    wfr.analysis(VERBOSE = True, DEBUG = True) 
    #print(wfr.custom_fields) ## good debug
        
    