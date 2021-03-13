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



###############################################################################
import sys
sys.path.append("/opt/WPG/") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/WPG") # DESY MAXWELL PATH

sys.path.append("../") # LOCAL PATH

###############################################################################
###############################################################################


from PIL import Image

import numpy as np
from matplotlib import pyplot as plt




def load_image_as_array(imdir, bPlot = False):
    """
    load a greyscale image as an array (for use as a phase or intenisty mask)
    
    :param imdir: directory of the input image
    
    :returns arr: array representation of surface height
    """
    
    arr = np.asarray(Image.open(imdir))
    
    if bPlot:
        plt.imshow(arr)
        
    return arr

def test():
    arr = load_image_as_array("../data/samples/AAO.png", bPlot = True)
    
    from wpg.srwl_uti_smp import srwl_opt_setup_transm_from_file as Sample
    
    #s = Sample(filepath, rx, ry, thickness, delta, atten_len, xc, yc, shift_x, shift_y)
    
    


if __name__ == '__main__':
    

    
    test()