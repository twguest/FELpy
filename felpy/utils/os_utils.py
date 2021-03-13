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


import os
import sys

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)
        
        
def np_size(nrows, ncols, ndepth, dtype=32, out="GB"):
   
   """Calculate the size of a numpy array in bytes.
   :param nrows: the number of rows of the matrix.
   :param ncols: the number of columns of the matrix.
   :param dtype: the size of each element in the matrix. Defaults to 32bits.
   :param out: the output unit. Defaults to gigabytes (GB)
   :returns: the size of the matrix in the given unit
   :rtype: a float
   """
   sizes = {v: i for i, v in enumerate("BYTES KB MB GB TB".split())}
   return nrows * ndepth * ncols * dtype / 8 / 1024. ** sizes[out] 


def add_path():
    """ 
    adds felpy filepath to python path
    """
    
    fpath = os.path.dirname(os.path.realpath(__file__))
    fpath = fpath.split("felpy/")[0]
    #print("Adding {} to Python path".format(os.path.join(fpath)))
    sys.path.append(fpath)
     
def felpy_path():
    """
    get felpy path
    """
    fpath = os.path.dirname(os.path.realpath(__file__))
    fpath = os.path.join(fpath.split("felpy/")[0] + "/felpy/")
    
    return fpath