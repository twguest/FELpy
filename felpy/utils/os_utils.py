#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "0.2.1"
__maintainer__ = "Trey Guest"
__email__ = "trey.guest@xfel.eu"
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
    if "/felpy/felpy/" in fpath:
        fpath = os.path.join(fpath.split("felpy/")[0] + "/felpy/felpy/")    
    else:
        fpath = os.path.join(fpath.split("felpy/")[0] + "/felpy/")

    return fpath

from functools import wraps
from time import time

def timing(f):
    @wraps(f)
    def wrap(*args, **kw):
        ts = time()
        result = f(*args, **kw)
        te = time()
        print('func:%r took: %2.4f sec' %(f.__name__,  te-ts))
        return result
    return wrap



class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = DEBUG = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = HEADER = '\033[4m'
   END = '\033[0m'

if __name__ == '__main__':
    print(felpy_path())
    