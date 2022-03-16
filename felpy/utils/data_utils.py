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



###############################################################################
import sys
sys.path.append("/opt/WPG/") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/WPG") # DESY MAXWELL PATH

sys.path.append("../") # LOCAL PATH

###############################################################################
###############################################################################

from os import listdir
from PIL import Image

import numpy as np
from matplotlib import pyplot as plt




def load_tif(imdir):
    """
    load a greyscale image as an array (for use as a phase or intenisty mask)
    
    :param imdir: location of the input image
    
    :returns arr: array representation of surface height
    """
        
    return np.asarray(Image.open(imdir))


def load_stack(imdir):
    """
    load a set of greyscale images as a 3 dimensional array
    
    :param imdir: directory containing input images
    
    :returns arr: x*y*n array of images
    """
    file_list = [file for file in listdir(imdir) if ".tif" in file]
    n = len(file_list)
    
    
    x,y = load_tif(imdir + listdir(imdir)[0]).shape
    
    stack = np.zeros([x,y,n])

    for itr in range(n):
        stack[:,:,itr] = load_tif(imdir + listdir(imdir)[itr])

    return stack


if __name__ == '__main__':
    load_stack("/media/twguest/Cache/imbl_data/sample_optimisation/proc/sic/1200grit_nocyl_1500ms/")