# -*- coding: utf-8 -*-
import numpy as np
import scipy
import scipy.ndimage as ndi
from matplotlib import pyplot as plt
def pshift(a, ctr):
    """
    Shift an array so that ctr becomes the origin.
    """
    sh  = np.array(a.shape)
    out = np.zeros_like(a)

    ctri = np.floor(ctr).astype(int)
    ctrx = np.empty((2, a.ndim))
    ctrx[1,:] = ctr - ctri     # second weight factor
    ctrx[0,:] = 1 - ctrx[1,:]  # first  weight factor

    # walk through all combinations of 0 and 1 on a length of a.ndim:
    #   0 is the shift with shift index floor(ctr[d]) for a dimension d
    #   1 the one for floor(ctr[d]) + 1
    comb_num = 2**a.ndim
    for comb_i in range(comb_num):
        comb = np.asarray(tuple(("{0:0" + str(a.ndim) + "b}").format(comb_i)), dtype=int)

        # add the weighted contribution for the shift corresponding to this combination
        cc = ctri + comb
        out += np.roll( np.roll(a, -cc[1], axis=1), -cc[0], axis=0) * ctrx[comb,range(a.ndim)].prod()

    return out

def free_nf(w, l, z, pixsize=1.):
    """\
    Free-space propagation (near field) of the wavefield of a distance z.
    l is the wavelength. 
    """
    if w.ndim != 2:
        raise RunTimeError("A 2-dimensional wave front 'w' was expected")

    sh = w.shape

    # Convert to pixel units.
    z = z / pixsize
    l = l / pixsize

    # Evaluate if alianp.sing could be a problem
    if min(sh)/np.sqrt(2.) < z*l:
        print("Warning: z > N/(sqrt(2)*lamda) = %.6g: this calculation could fail." % (min(sh)/(l*np.sqrt(2.)))) 
        print( "(consider padding your array, or try a far field method)"  )

    q2 = np.sum((np.fft.ifftshift(np.indices(sh).astype(float) - np.reshape(np.array(sh)//2,(len(sh),) + len(sh)*(1,)), range(1,len(sh)+1)) * np.array([1./sh[0], 1./sh[1]]).reshape((2,1,1)))**2, axis=0)

    return np.fft.ifftn(np.fft.fftn(w) * np.exp(2j * np.pi * (z / l) * (np.sqrt(1 - q2*l**2) - 1) ) )

  
# Simulation of a sphere
sh = (1024, 1024)
ssize = 5.    # rough speckle size
sphere_radius = 150
lam = .5e-10  # wavelength
z = 1.5    # propagation distance
psize = 1e-4  # pixel size

# Simulate speckle pattern
speckle = ndi.gaussian_filter(np.random.normal(size=sh), ssize) +\
          1j * ndi.gaussian_filter(np.random.normal(size=sh), ssize)



from .examples.NFS.cylinder_phase_mask import phase as cyl_phase

xx, yy = np.indices(sh)

xx, yy = np.indices(sh)
sphere = np.real(scipy.sqrt(sphere_radius**2 - (xx-256.)**2 - (yy-256.)**2))
sample = np.exp(-15*np.pi*2j*sphere/sphere_radius)

 
from felpy.experiments.speckle_tracking.optical_flow import process_optical_flow, kottler, process_all

Is = abs(free_nf(speckle*sample, 0.5e-10, z = .001, pixsize = psize))**2
Ir = abs(free_nf(speckle, 0.5e-10, z = .001, pixsize = psize))**2
plt.imshow(Is)
plt.show()
plt.imshow(Ir)
plt.show()

results = process_all(Is, Ir, sigma = .1, alpha = 1)

dx = results['dx']

dx[dx>1] = 1
dx[dx<-1]=-1


dy = results['dy']

dy[dy>1] = 1
dy[dy<-1]=-1

phi = results['phi']
from matplotlib import colors
plt.imshow(phi.real, vmax = 1)
#plt.imshow(dy.real)
