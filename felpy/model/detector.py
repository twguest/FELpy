# -*- coding: utf-8 -*-

""" 
Small detector module for scaling

"""

from felpy.model.mesh import Mesh
from felpy.model.beamline import Beamline
from wpg.optical_elements import Drift
from felpy.model.tools import propagation_parameters
import numpy as np
def rebin(arr, new_shape):
    """Rebin 2D array arr to shape new_shape by averaging."""
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)

def get_row_compressor(old_dimension, new_dimension):
    dim_compressor = np.zeros((new_dimension, old_dimension))
    bin_size = float(old_dimension) / new_dimension
    next_bin_break = bin_size
    which_row = 0
    which_column = 0
    while which_row < dim_compressor.shape[0] and which_column < dim_compressor.shape[1]:
        if round(next_bin_break - which_column, 10) >= 1:
            dim_compressor[which_row, which_column] = 1
            which_column += 1
        elif next_bin_break == which_column:

            which_row += 1
            next_bin_break += bin_size
        else:
            partial_credit = next_bin_break - which_column
            dim_compressor[which_row, which_column] = partial_credit
            which_row += 1
            dim_compressor[which_row, which_column] = 1 - partial_credit
            which_column += 1
            next_bin_break += bin_size
    dim_compressor /= bin_size
    return dim_compressor


def get_column_compressor(old_dimension, new_dimension):
    return get_row_compressor(old_dimension, new_dimension).transpose()

def compress_and_average(array, new_shape):
    # Note: new shape should be smaller in both dimensions than old shape
    return np.mat(get_row_compressor(array.shape[0], new_shape[0])) * \
           np.mat(array) * \
           np.mat(get_column_compressor(array.shape[1], new_shape[1]))
           
class Detector:
    
    def __init__(self, **kwargs):
        
        self.__dict__.update(kwargs)

        ### this means that we first defer to nx/ny/nz definitions of the field
        if all(hasattr(self, attr) for attr in ["nx", "ny", "dx", "dy"]):
            pass
        else:
            if hasattr(self, "mesh"):
                self.__dict__.update(self.mesh.get_attributes())
        
    def detect(self, wfr):
        """
        return the detected intensity distribution
        """
 
        
        bl = Beamline()
        bl.append(Drift(0), propagation_parameters(self.mesh.dx/wfr.dx, 1, self.mesh.dy/wfr.dy, 1))
        bl.propagate(wfr)
        
        ii = wfr.get_intensity().sum(-1)
        ii = compress_and_average(ii, (self.mesh.nx, self.mesh.ny))
        
        return ii
    

if __name__ == '__main__':
    
    from felpy.model.src.coherent import construct_SA1_pulse
    
    wfr = construct_SA1_pulse(512,512,3,10,0.25)
    [wxMin, wxMax, wyMax, wyMin] = wfr.get_limits()
    wnx, wny = wfr.nx, wfr.ny
    wdx, wdy = wfr.dx, wfr.dy
    

    
    detector_mesh = Mesh(nx = 400, ny = 270, dx = 5e-06, dy = 12e-06)
    
    d = Detector(mesh = detector_mesh)
    i = d.detect(wfr)
    
    print("Desired Properties/Actual Values")
    print("nx: {}/{}".format(detector_mesh.nx, wfr.params.Mesh.nx))
    print("ny: {}/{}".format(detector_mesh.ny, wfr.ny))
    print("dx: {}/{}".format(detector_mesh.dx, wfr.dx))
    print("dy: {}/{}".format(detector_mesh.dy, wfr.dy))