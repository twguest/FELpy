# -*- coding: utf-8 -*-


from wpg.srwl_uti_smp import srwl_opt_setup_transm_from_file as Sample
from felpy.utils.vis_utils import colorbar_plot
from PIL import Image
from felpy.utils.np_utils import get_mesh
import matplotlib.colors as colors
from matplotlib import pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from felpy.utils.opt_utils import ekev2wav, get_phase_shift, wrap_phase
import numpy as np


def rgb2gray(rgb):

    r, g, b = rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]
    gray = 0.2989 * r + 0.5870 * g + 0.1140 * b

    return gray

def load_tif(directory):
    """
    load a tif file as an np array
    """
    return (np.asarray(Image.open(directory)))
    
def define_speckle_mask(sample_directory, rx, ry, sample_thickness, sample_delta, ekev = None, xc = 0, yc = 0, plot = False):
    
    arr = load_tif(sample_directory)
    
    sample_thicc = arr/255*sample_thickness

    if plot:
        plt.style.use(['science','ieee'])
        
        
        ### plot sample thickness
        
 
        fig, ax1 = plt.subplots()
        img = ax1.imshow(sample_thicc*1e6, cmap = "bone",norm=colors.Normalize())

        ax1.set_xlabel("x (mm)")
        ax1.set_ylabel("y (mm)")

        divider = make_axes_locatable(ax1)
        cax = divider.append_axes('right', size='7.5%', pad=0.05)
        
        cbar = fig.colorbar(img, cax)
        cbar.set_label("Thickness ($\mu$m)")
    
    if plot and ekev is not None:
        plt.style.use(['science','ieee'])
        
        
        ### plot sample phase-shift
        

        fig, ax1 = plt.subplots()
        img = ax1.imshow(wrap_phase(get_phase_shift(sample_thicc,ekev2wav(ekev))), cmap = "hsv")
        

        ax1.set_xlabel("x (mm)")
        ax1.set_ylabel("y (mm)")

        divider = make_axes_locatable(ax1)
        cax = divider.append_axes('right', size='7.5%', pad=0.05)
        
        cbar = fig.colorbar(img, cax)
        cbar.set_label("Phase Shift (rad)")
        
        cax.set_xlim([-np.pi, np.pi])
        
        
    s = Sample(sample_directory, 
           rx = rx, ry = ry,
           thickness = sample_thicc,
           delta = sample_delta,
           atten_len = 100,
           xc = xc, yc = yc,
           invert = True)
    
    
    return s


if __name__ == '__main__':
    import numpy as np
    arr = load_tif("../../data/samples/speckle_enlarged.tif")
    
    
    sample = define_speckle_mask("../../data/samples/speckle.tif", 1, 1, sample_thickness = 21e-06,
                                 sample_delta = 4e-07, plot = True, ekev = 25)