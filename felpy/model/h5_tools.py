#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "1.0.1"
__maintainer__ = "Trey Guest"
__email__ = "twguest@students.latrobe.edu.au"
__status__ = "Developement"
"""

import h5py
import numpy as np
import datetime

class obj:
    
    def __init__(self):
        self.created = str(datetime.datetime.now())
        self.author = "twguest"
        pass


def object2h5(fname, obj):
    """ 
    write a python object to hdf5 file
    
    :param fname: filename for h5 writing
    :param obj: object to be written
    """
    with h5py.File(fname, 'w') as f:
        for item in vars(obj).items():
            
            f.create_dataset(item[0], data = item[1])
    

def h52object(fname):
    """
    read a hdf5 file to a python object
    
    :param fname: filename for h5 writing
    """
    o = obj()

    with h5py.File(fname, 'r') as f:
        for key in f.keys():
            setattr(o, key, f[key].value)
            
            
    return o



if __name__ == '__main__':
    
    
    ### test
    
    pyobj = obj()
    
    ### create some random object components
    
    pyobj.a = "a"
    pyobj.b = [1,2,3,4,5]
    pyobj.c = np.random.rand(5,5)
    pyobj.d = np.random.rand(5,5) #+ np.random.rand(5,5)*1j
    
    ### test object to h5
    
    object2h5("testObject.h5", pyobj)
    
    ### test h5 to object
    
    obj = h52object("testObject.h5")
