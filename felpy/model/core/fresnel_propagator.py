# manage the imports
import numpy as np
 
 
def fft(x):
    # this only defines the correct fwd fourier transform including proper shift of the frequencies
    return np.fft.fftshift(np.fft.fft2(x)) # Fourier transform and shift

def ifft(x):
    return np.fft.ifft2(np.fft.ifftshift(x)) # inverse Fourier transform and shift
    
def frensel_propagator(wfr, px, py, wav, z, upsample = 1):
    """
    fresnel propagator
    
    :param wfr: complex wfr [nx,ny] 
    :param px: horizontal pixel size [m]
    :param py: vertical pixel size [m]
    :param wav: wavelength [m]
    :param z: distance to propagate [m]
    :param upsample: upsampling factor of grid
    
    :returns Ef: propagated field
    """

    grid_size_x = px * upsample;                 # Grid size in x-direction
    grid_size_y = py * upsample;                 # Grid size in x-direction
      

    # Inverse space
    fx = np.linspace(-(upsample-1)/2*(1/grid_size_x), (upsample-1)/2*(1/grid_size_x), upsample)
    fy = np.linspace(-(upsample-1)/2*(1/grid_size_y), (upsample-1)/2*(1/grid_size_y), upsample)
    Fx, Fy = np.meshgrid(fx, fy)
    
 
    H = np.exp(1j*(2 * np.pi / wav) * z) * np.exp(1j * np.pi * wav * z * (Fx**2 + Fy**2))
     
    # Compute FFT centered about 0
    wfrfft = fft((wfr));     # Centered about 0 since fx and fy centered about 0
    
    # Multiply spectrum with fresnel phase-factor
    G = H * wfrfft
    Ef = ifft(G) # Output after deshifting Fourier transform
    
    return np.moveaxis(Ef, -1, 0)
