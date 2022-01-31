#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "0.1.1"
__maintainer__ = "Trey Guest"
__email__ = "twguest@students.latrobe.edu.au"
__status__ = "Developement"
"""

import h5py
import types
import numpy as np


def blank_object():
    """
    construct a blank python object
    """
      
    obj = types.SimpleNamespace()
    
    return obj

    
def create_h5(sdir):
    """ 
    create a h5 file
    
    :param sdir: save directory (incl. filename) of .h5 file
    """    
    with h5py.File(sdir + ".hdf5", 'w') as f:
        f.close()
        
 
        
def create_cxi_format(sdir):
    """
    create a .h5 file in the CXI format
    
    :param sdir: save directory (incl. filename) of .h5 file
    """
    
    create_h5(sdir)
    
    with h5py.File(sdir + ".hdf5", 'w') as f:
        
        f.create_group("/entry_1/")
        
        f.create_group("/entry_1/data_1/")
        
        f.create_dataset("/entry_1/data_1/data",
                         data = [])

        f.create_group("/entry_1/instrument_1/")
        
        f.create_group("/entry_1/instrument_1/detector_1/")
        
        f.create_dataset("/entry_1/instrument_1/detector_1/x_pixel_size",
                         data = [])      
        
        f.create_dataset("/entry_1/instrument_1/detector_1/y_pixel_size",
                         data = [])   
        
        f.create_dataset("/entry_1/instrument_1/detector_1/basis_vectors",
                         data = [])   
        
        f.create_dataset("/entry_1/instrument_1/detector_1/mask",
                         data = [])   
        
        f.create_dataset("/entry_1/instrument_1/detector_1/distance",
                         data = [])        
       
        f.create_group("/entry_1/source_1/")
        
        f.create_dataset("/entry_1/instrument_1/source_1/energy",
                         data = [])        
        
        f.create_dataset("/entry_1/instrument_1/source_1/wavelength",
                         data = [])   
        
        f.create_group("/entry_1/sample_1/")
        
        f.create_group("/entry_1/sample_1/geometry/")
        
        f.create_dataset("/entry_1/sample_1/geometry/translation",
                         data = [])      
    
        f.close()
        
        
if __name__ == '__main__':
    
    sdir = "/tmp/practice_file"
    
    create_cxi_format(sdir)