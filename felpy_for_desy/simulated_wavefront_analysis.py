# -*- coding: utf-8 -*-


import sys

from felpy.model.core.wavefront import Wavefront

if __name__ == '__main__':
    
    in_directory = sys.argv[1]
    out_directory = sys.argv[2]
    
    wfr = Wavefront()
    wfr.analysis()
    
        
    