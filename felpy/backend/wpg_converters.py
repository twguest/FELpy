# -*- coding: utf-8 -*-

import numpy as np

from felpy.model.wavefront import Wavefront


def complex_to_wpg(arr): ### converter
    """
    converter function to transform complex wavefield into wpg style electric
    field array
    
    :param arr: complex wavefield array [x,y,t] (complex128 type)
    
    :returns new_arr: wpg style electric field array [nx,ny,nz,2] (float64)
    """
    new_arr = np.zeros([arr.shape[0], arr.shape[1], arr.shape[2], 2])
    new_arr[:,:,:,0] = arr.real
    new_arr[:,:,:,1] = arr.imag
    return new_arr



def wavefront_from_array(cfr,nx,ny,nz,dx,dy,dz,ekev, pulse_duration = 40e-15, sigma = 4):
    """
    function to produce a wpg wavefront object instance from a complex valued
    wavefront definition
    
    :param cfr: complex valued wavefield array [x,y,t] (complex128)
    :param nx: number of pixels along x-axis
    :param ny: number of pixels along y-axis
    :param nz: number of pixels along z-axis
    :param dx: pixels size along x-axis
    :param dy: pixels size along y-axis
    :param dz: pixels size along z-axis
    :param ekev: energy in kev
    :param pulse_duration: duration of the pulse (fwhm) in seconds, defaults to 40e-15
    :param sigma: integer number of width multiples over which time-axis should be defined, defaults to 4

    :returns: DESCRIPTION
    """
    
    # Initialize empty wavefront.
    wfr = Wavefront()

    # Setup E-field.
    wfr.data.arrEhor = np.zeros(shape=(nx, ny, nz, 2))
    wfr.data.arrEver = np.zeros(shape=(nx, ny, nz, 2))

    wfr.params.wEFieldUnit = 'sqrt(W/mm^2)'
    wfr.params.photonEnergy = ekev * 1000
    wfr.params.wDomain = 'time'
    wfr.params.Mesh.nSlices = nz
    wfr.params.Mesh.nx = nx
    wfr.params.Mesh.ny = ny


    
    wfr.params.Mesh.sliceMin = -pulse_duration*sigma / 2.
    wfr.params.Mesh.sliceMax = pulse_duration*sigma / 2.

    range_x = dx*nx
    range_y = dy*ny

    wfr.params.Mesh.xMin = -range_x / 2.
    wfr.params.Mesh.xMax = range_x / 2.
    wfr.params.Mesh.yMin = -range_y / 2.
    wfr.params.Mesh.yMax = range_y / 2.

    wfr.params.Rx = 2
    wfr.params.Ry = 1
    
    wfr.data.arrEhor = complex_to_wpg(cfr)
    
    wfr.set_electric_field_representation('f')
        
    return wfr