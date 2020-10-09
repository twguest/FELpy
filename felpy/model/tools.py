
###############################################################################
import sys
sys.path.append("/opt/WPG/") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/WPG") # DESY MAXWELL PATH

sys.path.append("/opt/spb_model") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/spb_model") # DESY MAXWELL PATH
###############################################################################
###############################################################################
import numpy as np
import os
from wpg import srwlib
from wpg.wavefront import Wavefront
from wpg.wpg_uti_wf import integral_intensity, calculate_fwhm, look_at_q_space
from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity
from felpy.model.src.coherent import coherentSource
from wpg.srwlib import SRWLOptD as Drift
from wpg.srwlib import SRWLOptA as Aperture
from wpg.generators import build_gauss_wavefront
from wpg.beamline import Beamline

def propParams(sx, zx, sy, zy, mode = "normal"):
    """
    wrapper for propagation parameters
    
    :param sx: horizontal scaling factor 
    :param zx: horizontal zoom factor
    :param sy: vertical scaling factor
    :param zy: vertical zoom factor
    :param mode: normal, semi-analytical, converge or diverge
    
    :return propagation parameters:
    """
    
    if mode == "fresnel":
        m = 0
    elif mode == "quadratic":
        m = 1
    elif mode == "fraunhofer":
        m = 2
    elif mode == "diverge":
        m = 3
    elif mode == "converge":
        m = 4
    
    return [0,0,1,m,0,sx,zx/sx,sy,zy/sy,0,0,0]

def constructPulse(nx = 512, ny = 512, nz = 5, tau = 1e-06, d2waist = 10):
    
    wfr = Wavefront(build_gauss_wavefront(nx, ny, nz, 5.0, -400e-06, 400e-06, -400e-06, 400e-06, tau, 5e-06, 5e-06, d2waist))
    srwlib.srwl.SetRepresElecField(wfr._srwl_wf, 'f')
    #look_at_q_space(wfr)
    return wfr

def create_circular_mask(nx, ny, c = None, r = None):
    """
    create a circular mask numpy array
    
    :param nx: number of pixels of array (horizontal) [int]
    :param ny: number of pixels of array (vertical) [int]
    :param c: center of array in pixels [list]from wpg import srwlib
    :param r: radius of array (in pixels) [int]
    
    :returns mask: binary mask [np array]
    """ 
    if c is None: # use the middle of the image
        c = (int(nx/2), int(ny/2))
    if r is None: # use the smallest distance between the center and image walls
        r = min(c[0], c[1], nx-c[0], ny-c[1])

    X, Y = np.ogrid[:nx, :ny]
    dist_from_center = np.sqrt((X - c[0])**2 + (Y-c[1])**2)

    mask = dist_from_center <= r
    
    return mask

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)
        

def radial_profile(data, center):
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile, r


def binArray(data, axis, binstep, binsize, func=np.nanmean):
    data = np.array(data)
    dims = np.array(data.shape)
    argdims = np.arange(data.ndim)
    argdims[0], argdims[axis]= argdims[axis], argdims[0]
    data = data.transpose(argdims)
    data = [func(np.take(data,np.arange(int(i*binstep),int(i*binstep+binsize)),0),0) for i in np.arange(dims[axis]//binstep)]
    data = np.array(data).transpose(argdims)
    return data


def getCoord(arr, mode = 'centre'):
    
    """
    get the coordinate space indices of a 2D array, 2 2D array of indices.
    useful in translating numpy coordinate system to real-space coordinates
    
    :param arr: np array to get index of
    :param mode: if center, return indices relative to the matrix center (0,0) = matrix center,
                 otherwise, return relative to numpy representation (0,0) = top left corner 
    
    :usage:
        >>> getIndex(arr = np.ones([9,9]), mode = 'centre')
        
    :returns idx: x-coordinates (2d numpy array)
    :returns idy: y-coordinates (2d numpy array)
    """
    
    if (arr.shape[0]%2) == 0 or (arr.shape[1]%2) == 0:
        raise(ValueError("Number of Rpeates Should Be Odd"))
    
    
    if mode == 'centre':
        idx = np.linspace(0,arr.shape[0], arr.shape[0]+1)[1:]-(arr.shape[0]//2+1)
        idy = np.flip(np.linspace(0,arr.shape[1], arr.shape[1]+1)[1:]-(arr.shape[1]//2+1))
    else:
        idx = np.linspace(0,arr.shape[0], arr.shape[0]+1)[1:]
        idy = np.linspace(0,arr.shape[1], arr.shape[1]+1)[1:]
        
    idx, idy = np.meshgrid(idx, idy)
    
    return idx, idy 

def argmax2d(X):
    """
    2-dimensional equivalent of np.argmax
    
    via: https://github.com/numpy/numpy/issues/9283
    """
    n, m = X.shape
    x_ = np.ravel(X)
    k = np.argmax(x_)
    i, j = k // m, k % m
    return i, j

def memoryMap(savedir, fname, shape, dtype = 'float64'):
    """
    construct a memory map
    """
    if os.path.exists(savedir + fname):
        memmap = np.memmap(savedir + fname, mode = "r+",
                           shape = shape, dtype = dtype)
    else:
        memmap = np.memmap(savedir + fname, mode = "w+",
                           shape = shape, dtype = dtype)
    return memmap

def readMap(mapdir,shape,dtype = 'float64'):
    """
    read a map from mapDir
    """
    
    mp = np.memmap(mapdir, dtype = dtype, mode = 'r+', shape = shape)
    
    return mp


def generateTestPulses(savedir, nx = 1024, ny = 1024, N = 5):
    """
    generate a set of test pulses
    
    :param savedir: directory for test pulses to be saved
    :param N: number of pulses to generate
    """
    print("Constructing {} Test Pulses".format(N))
    
    mkdir_p(savedir)
    
    for n in range(N):
        
        wfr = constructPulse(nx, ny, nz = 6, tau = 1e-12)
        
        wfr.data.arrEhor*= np.random.uniform(0.75, 1, size = wfr.data.arrEhor.shape)
        
        wfr.store_hdf5(savedir + "testWavefront_{}.h5".format(n))
        print("Storing Wavefront @" + savedir + "testWavefront_{}.h5".format(n))    