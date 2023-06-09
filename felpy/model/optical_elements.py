import numpy as np
from wpg.srwlib import srwl_opt_setup_surf_height_2d as MirPl
from wpg.optical_elements import Mirror_plane_2d

from scipy import interpolate
from wpg.srwlib import SRWLOptT
import os

def Mirror_plane_2d(orient, theta, length, range_xy, _height_prof_data, xpos, ypos, refl = 1,  scale=1, x0=0., y0=0., xscale=1., yscale=1.):
    """
    Defining a plane mirror propagator with taking into account 2D surface height errors

    :param orient:  mirror orientation, 'x' (horizontal) or 'y' (vertical)
    :param theta:   incidence angle [rad]
    :param length:  mirror length, [m]
    :param range_xy: range in which the incident WF defined [m]
    :param _height_prof_data:2d mirror profile h(x, y) - heigh errors [m]
    :param x:
    :param y:
    :param refl:
    :param scale: scale factor, optical path difference OPD = 2*h*scale*sin(theta)
    :param x0: shift of mirror longitudinal position [m]
    :param y0: shift of mirror transverse position [m]
    :param xscale: units of 1st column of filename,  x[m]=x[nits]*xscale  [m]
    :param yscale: units of 1st column of filename,  y[m]=y[nits]*yscale  [m]
    :return: opIPM  - imperfect plane mirror propagator
    """


    sinTheta = np.sin(theta)
 
    
    dim = np.shape(_height_prof_data)
    ntotal = dim[0]
    
    nx, ny = _height_prof_data.shape
    print('nx,ny:', nx, ny)

 
    xmin = min(xpos)
    xmax = max(xpos)
    xc = (xmin+xmax)/2
 
    ymin = min(ypos)
    ymax = max(ypos)
    yc = (ymin+ymax)/2

    xpos = xpos - x0 - xc
    
    xmin = min(xpos)
    xmax = max(xpos)

    ypos = ypos - y0 - yc
    ymin = min(ypos)
    ymax = max(ypos)

    print('length: {:.1f} mm, width: {:.1f} mm'.format(
        (xmax-xmin)*1e3, (ymax-ymin)*1e3))

    # if (xmin <= -length/2.) and (xmax >= length/2):
    #     xmin = -length/2
    #     xmax = length/2
    # else:
    #     raise ValueError(
    #         'specified length exceeds mirror limits')
    # if (ymin <= -range_xy/2) and (ymax >= range_xy/2):
    #     ymin = -range_xy/2
    #     ymax = range_xy/2
    # else:
    #     raise ValueError(
    #         'specified width exceeds mirror limits'
    #     )


    print(xmax, xmin)
    print(ymax, ymin)
    if orient == 'y':
        opIPM = SRWLOptT(nx, ny, (ymax-ymin), (xmax-xmin)*sinTheta,_x = xc, _y = yc)
    elif orient == 'x':
        opIPM = SRWLOptT(nx, ny, (xmax-xmin)*sinTheta, (ymax-ymin),_x = xc, _y = yc)
    else:
        raise TypeError('orient should be "x" or "y"')
    xnew, ynew = np.mgrid[xmin:xmax:nx*1j, ymin:ymax:ny*1j]
    f = interpolate.RectBivariateSpline(xpos, ypos, _height_prof_data.T)
    h_new = f(xnew[:, 0], ynew[0, :])
 
    auxMesh = opIPM.mesh
    
    from array import array
    foo = array(str(u'd'), [])
    # for i in range(150000):
    #     foo.append(1.)

    foo = array(str(u'd'), [1.]*nx*ny)
    
    opIPM.arTr[::2] = foo  # Amplitude Transmission
 
    foo = array(str(u'd'), [])
    
    if orient == 'y':
        for ix in range(nx):
            for iy in range(ny):
                foo.append(-2 * sinTheta * h_new[ix, iy] * scale)
   
    elif orient == 'x':
        for iy in range(ny):
            for ix in range(nx):
                foo.append(-2 * sinTheta * h_new[ix, iy] * scale)

    opIPM.arTr[1::2] = foo    # Optical Path Difference (to check sign!)

    return opIPM