# -*- coding: utf-8 -*-

from felpy.model.src.coherent import construct_SA1_wavefront
import numpy as np
from felpy.utils.np_utils import get_mesh, equate_mesh, crop_to
from felpy.utils.vis_utils import basic_plot

def generate_wfr(nx, ny):
    
    return construct_SA1_wavefront(nx, ny, 5.0, 0.25)


def generate_exp(nx, ny):
    
    return np.random.rand(nx, ny)
    
def compare_simulation(wfr, ii, mesh):
    """
    A function to compare the intensity distributions of a wpg.wavefront
    object and an experimental wavefront. Crops the larger of the two-objects
    to the size of the smaller.
    
    :param wfr: simulation wavefront (wpg.wavefront)
    :param ii: experimental intensity distribution
    :param mesh: experimental mesh (ie get_mesh(ii, dx, dy))
    """
    
    wpg_ii = wfr.get_intensity().sum(-1)
    wpg_mesh = wfr.get_mesh()
    
    if np.max(wpg_mesh) > np.max(mesh):
        wpg_mesh = equate_mesh(mesh, wpg_mesh)
    elif np.max(mesh) > np.max(wpg_mesh):
        mesh = equate_mesh(wpg_mesh, mesh)
    
    ii, wpg_ii = crop_to(ii, wpg_ii)
    
    return [wpg_ii, wpg_mesh], [ii, mesh]
        

if __name__ == '__main__':
    
    nx = ny = 512
    px = py = np.random.randint(5, 150)
    
    efr = generate_exp(nx-px, ny-py)    
    sfr = generate_wfr(nx, ny)
    sx, sy = sfr.get_spatial_resolution()
    emesh = get_mesh(efr, sx, sy)
    exp, sim = compare_simulation(sfr, efr, emesh)
    
    basic_plot(*sim, cmap = "jet")
    basic_plot(*exp, cmap = "jet")