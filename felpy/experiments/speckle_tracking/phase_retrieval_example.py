import os

import numpy as np
from felpy.utils.np_utils import load_tif

from felpy.experiments.speckle_tracking.paganin_method import paganin_method

from felpy.experiments.speckle_tracking.optical_flow import process_optical_flow
from felpy.experiments.speckle_tracking.optical_flow import kottler as kottler
from felpy.experiments.speckle_tracking.optical_flow import LarkinAnissonSheppard as las
from felpy.experiments.speckle_tracking.frankoChellappa import frankotchellappa as franko

from felpy.experiments.speckle_tracking.fast_phase_retrieval import geometric_flow

from felpy.experiments.speckle_tracking.UMPA import match_speckles as UMPA
from felpy.experiments.speckle_tracking.UMPA import get_measurements

from felpy.utils.opt_utils import ekev2wav
from felpy.utils.vis_utils import Grids

from felpy.utils.np_utils import get_mesh, extent_from_mesh


def plot_images(reference_map, distorted_maps, cmap, clabel):
    
    mesh = get_mesh(reference_map, 1, 1)
    
    extent = [np.min(mesh[1]), np.max(mesh[1]),
              np.min(mesh[0]), np.max(mesh[0])]
    
    n = 1
    m = 4
    
    vmin = np.min(reference_map.real)
    vmax = np.max(reference_map.real)
    
    plots = Grids(global_aspect = 5, scale = 2)
    plots.create_grid(n = 1, m = 4, xlabel = "x (mm)",
                                    ylabel = "y (mm)",
                                    title = "")

    
    plots.axes[0].imshow(reference_map.real, cmap = cmap, extent = extent)
    
    for i in range(3):
        plots.axes[i+1].imshow(distorted_maps[i].real,
                                                    cmap = cmap,
                                                    extent = extent)
        
    plots.add_global_colorbar(clabel = clabel, cmap = cmap,
                              vmin = int(vmin), vmax = int(vmax))
    plots.fig.show()


def determine_thickness(phase_shift, wav, delta):
    """ 
    determine the projected thickness of a material from its prospective phase
    shift
    """  
    return (wav*phase_shift)/(2*np.pi*delta)


if __name__ == '__main__':
    ### sample image loading/generation
    
# =============================================================================
#     
#     ### using sample data
#     fdir = "../../data/zyla_speckle/" ### file directory
#     images = [load_tif(fdir + a)[700:800,700:800] for a in os.listdir(fdir)]
#     #images = [load_tif(fdir + a) for a in os.listdir(fdir)]
#     
#     reference = images[0]
#     images = images[1:]
# =============================================================================
     
    ### using thibault speckle
    reference, images = get_measurements()
    images = list(images)
    
    
    plot_images(reference, images, cmap = 'bone', clabel = "Projected Thickness")

    
    reference_map = paganin_method(I = reference, n_ratio = 100,
                                   wav = ekev2wav(6.2), z = .15)
    
    distorted_maps = [paganin_method(I = a, n_ratio = 100, wav = ekev2wav(6.2),
                                     z = .15)
                      for a in images]
    
    reference_thickness = determine_thickness(reference_map, ekev2wav(6.2), delta = 4e-07)
    
    distorted_thickness = []
    for item in distorted_maps:
        distorted_thickness.append(determine_thickness(item, ekev2wav(6.2), delta = 4e-07))
    
    
    plot_images(reference_thickness, distorted_thickness, cmap = 'viridis', clabel = "Projected Thickness")
    
    displacements = []
    reference_shift = np.zeros_like(reference_map)  
    
    for item in distorted_maps:
        displacements.append(process_optical_flow(item.real, reference_map.real))
    
    dx = [displacements[i][0] for i in range(len(displacements))]
    dy = [displacements[i][1] for i in range(len(displacements))]
    
    plot_images(reference_shift, dx, cmap = 'coolwarm', clabel = "Hor. Displacement (pix.)")
    plot_images(reference_shift, dy, cmap = 'coolwarm', clabel = "Ver. Displacement (pix.)")
    
    reference_phase_shift = fast_phase_retrieval(reference_shift, reference_shift, ekev2wav(6.2), z = .15)
    
    phase_shifts = []
    for itr in range(len(dx)):
        phase_shifts.append(fast_phase_retrieval(dx[itr], dy[itr], wav = ekev2wav(6.2), z = .15))
    
    plot_images(reference_shift, phase_shifts, cmap = 'bone', clabel = "Projected Phase Shift")
