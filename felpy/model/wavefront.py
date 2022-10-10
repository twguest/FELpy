# -*- coding: utf-8 -*-
"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "0.2.1"
__maintainer__ = "Trey Guest"
__email__ = "trey.guest@xfel.eu"
__status__ = "Developement"
"""


from wpg.srw import srwlpy 
from wpg.wavefront import Wavefront as WPG_Wavefront
from matplotlib import pyplot as plt
import imageio
import numpy as np
from wpg.wpg_uti_wf import calculate_fwhm

import numpy as np
from wpg.wpg_uti_wf import calculate_fwhm, calc_pulse_energy
import seaborn as sns
from felpy.utils.np_utils import get_mesh
from felpy.analysis.enclosed_energy import get_enclosed_energy
from scipy.constants import h, c, e
from felpy.analysis.centroid import get_com
from felpy.analysis.complex.coherence import get_coherence_time


from felpy.utils.np_utils import get_mesh
from felpy.model.tools import radial_profile
from datetime import datetime
import os
from felpy.analysis.energy_spectral_density import power_spectral_density
from felpy.analysis.centroid import get_com
from felpy.utils.maths.fit import fit_gaussian
from felpy.utils.maths.constants import sigma_to_fwhm

ls = {"m": 1,
      "cm": 1e2,
      "mm": 1e3,
      "um": 1e6,
      "nm": 1e9}

def complex_converter(carr):
    """
    convert a complex array to wpg format (note, currently only taking hor.
                                           polarised component).
    
    :param carr: numpy style complex array
    """
    if len(carr.shape) == 3:
        
        cwfr = np.ones([carr.shape[0], carr.shape[1], 2, carr.shape[2]])
        cwfr[:,:,0,:] = carr.real
        cwfr[:,:,1,:] = carr.imag
    else:
        
        cwfr = np.ones([carr.shape[0], carr.shape[1], 1,2])
    
        cwfr[:,:,0,0] = carr.real
        cwfr[:,:,0,1] = carr.imag

    return cwfr


class Wavefront(WPG_Wavefront):
    
    def __init__(self, _srwl_wf = None):
        super().__init__(_srwl_wf)
        self.custom_fields['metadata'] = {}
        self._changed = False

            
    def load_complex_array(self, carr):
        """
        load a complex array in WPG format, assumes no change in f.o.v etc.
        """
        self.data.arrEhor /= self.data.arrEhor 
        self.data.arrEhor *= complex_converter(carr)
        
    def multiply_by_complex(self, carr):
        """
        multiply the current electric field by a complex array
        """
        self.data.arrEhor[:,:,:,0] *= (complex_converter(carr)[:,:,:,0]).astype(np.float32)
        self.data.arrEhor[:,:,:,1] += (complex_converter(carr)[:,:,:,1]).astype(np.float32)
        
    def save_tif(self, outdir):
        ii = self.get_intensity().sum(-1)
        imageio.imwrite(outdir + ".tif", ii)
 
    def get_spatial_resolution(self, VERBOSE = False):
        
        px = (self.params.Mesh.xMax - self.params.Mesh.xMin) / self.params.Mesh.nx
        py = (self.params.Mesh.yMax - self.params.Mesh.yMin) / self.params.Mesh.ny
        
        
        self.custom_fields['spatial resolution'] = px, py
        
        
        if VERBOSE:
            print(self.custom_fields['spatial resolution'])
        return px, py
    
    def get_angular_resolution(self, VERBOSE = False):
        qx = (self.params.Mesh.qxMax - self.params.Mesh.qxMin) / self.params.Mesh.nx
        qy = (self.params.Mesh.qyMax - self.params.Mesh.qyMin) / self.params.Mesh.ny
        
        
        self.custom_fields['angular resolution'] = qx, qy
        
        
        if VERBOSE:
            print(self.custom_fields['angular resolution'])
        return qx, qy
    
    @property
    def pulse_duration(self):
        
        wDomain = self.params.wDomain
        
        if wDomain == 'time':
            pulse_duration = (self.params.Mesh.sliceMax - self.params.Mesh.sliceMin) 
        else:
            self.set_electric_field_representation('t')
            pulse_duration = (self.params.Mesh.sliceMax - self.params.Mesh.sliceMin) 
            self.set_electric_field_representation(wDomain)
            
 
        return pulse_duration


    def to_txt(self, outdir):
        
  
        output = open(outdir + ".txt", "w")
        output.write(self.__str__())
        output.close()
    
    
    
    @property
    def spatial_power_spectral_density(self):
        """
        wrapper function for felpy.analysis.energy_spectral_density.power_spectral_density
        """
        freq, s = power_spectral_density(self.as_complex_array().sum(-1),
                                         spatial_sampling = self.dx, pulse_duration = self.pulse_duration)
        
        return freq, s
    def as_complex_array(self, ignoreVer = True, ignoreHor = False):
        """
        Convert electric field data to complex representation]
        
        returns [Ehor, Ever]
        """
        
        
        for slc in range(self.params.Mesh.nSlices):
            
                       
            if slc == 0:
                
                Ehor = self.data.arrEhor[:,:,slc,0].astype('complex128')+(self.data.arrEhor[:,:,slc,1]*1j)
                Ehor = Ehor.reshape((*self.data.arrEhor.shape[0:2], 1))
                
                Ever = self.data.arrEver[:,:,slc,0].astype('complex128')+(self.data.arrEver[:,:,slc,1]*1j)
                Ever = Ever.reshape((*self.data.arrEver.shape[0:2], 1))
            else:
                Etmp = (self.data.arrEhor[:,:,slc,0].astype('complex128')+(self.data.arrEhor[:,:,slc,1]*1j)).reshape((*self.data.arrEhor.shape[0:2], 1))
                Ehor = np.concatenate((Ehor, Etmp), axis = 2)
                
                Etmp = (self.data.arrEver[:,:,slc,0].astype('complex128')+(self.data.arrEver[:,:,slc,1]*1j)).reshape((*self.data.arrEver.shape[0:2], 1))
                Ever = np.concatenate((Ever, Etmp), axis = 2)
                
                
        Ehor.imag, Ever.imag = np.imag(Ehor), np.imag(Ever)
        
        return np.array([Ehor, Ever])[0,:,:,:]### PRETTY SURE THIS JUST RETURNS Ehor
    
    
    def view(self):
        
        """
        useless - should be deprecated, just use wpg utils
        
        a simple method for viewing the intensity and phase of a wavefront
        if a deeper analysis is required. try wpg.wpg_uti_wf.plot_intensity_map()
        """
        fig, axs = plt.subplots(1,2)
        [ax1, ax2] = axs
        
        ax1.imshow(self.get_intensity().sum(axis = -1), cmap = 'bone')
        ax1.set_title("Mean Intensity")
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax2.set_xticks([])
        ax2.set_yticks([])
        ax2.imshow(self.get_phase().sum(axis = -1), cmap = 'hsv')
        ax2.set_title("Mean Phase")
        plt.show()


    def get_wavelength(self, VERBOSE = False):
        
        self.custom_fields['wavelength'] = (h*c)/(e*self.params.photonEnergy)
        
        if VERBOSE: print("wavelength")
        if VERBOSE: self.custom_fields['temporal resolution']
        
        return self.custom_fields['wavelength'] 


    def set_electric_field_representation(self, domain):
        """
        wrapper for srwlpy.SetRepresElecField
        
        sets the electric field representation
        
        :param domain: choice ofangular - 'a', frequency 'f' or time 't' [str]
        """
        srwlpy.SetRepresElecField(self._srwl_wf, domain)
        

    def get_divergence(self):
        """
        calculate the full-angle divergence of the beam
        
        :param wfr: WPG wavefront structure
        """
        
        self.set_electric_field_representation('a') 
        sig_x, sig_y = calculate_fwhm(self)['fwhm_x'], calculate_fwhm(self)['fwhm_y']
        
        self.set_electric_field_representation('c') 
        
        self.custom_fields['divergence'] = sig_x, sig_y
        
        return sig_x, sig_y
    
    def get_fwhm(self):
        sig_x, sig_y = calculate_fwhm(self)['fwhm_x'], calculate_fwhm(self)['fwhm_y']
        
        self.custom_fields['fwhm'] = sig_x, sig_y
        
        return sig_x, sig_y


    def get_energy_statistics(self, integrate = False, VERBOSE = False, mpi = True, write = False):
        

        self.set_electric_field_representation('t')
        energy = calc_pulse_energy(self)
        self.set_electric_field_representation('f')
    
        self.custom_fields['pulse energy'] = energy[0]
        self.custom_fields['nphotons'] = energy[1]
        
        return energy 
     
        
    def get_pulse_duration(self, VERBOSE = False):
        
        self.set_electric_field_representation('t')
        t = self.params.Mesh.nSlices*self.get_temporal_resolution()
        self.set_electric_field_representation('f')
        
        self.custom_fields['pulse duration'] = t
        
        if VERBOSE: 
            print("Pulse Duration")
            print(self.custom_fields['pulse duration'])
            
        return t
    
    def plot_intensity(self, scale = "mm", label = "", title = "", context = 'talk', sdir = None):
        
        ii = self.get_intensity().sum(-1)
        
        ls = {"m": 1,
              "cm": 1e2,
              "mm": 1e3,
              "um": 1e6,
              "nm": 1e9}
        
        mesh = get_mesh(ii, *self.get_spatial_resolution())

        colorbar_plot(ii, mesh,
                      label = label,
                      title = title,
                      xlabel = "x ({})".format(scale),
                      ylabel = "y ({})".format(scale),
                      clabel = "Intensity (a.u.)",
                      vmax = np.max(ii),
                      context = context,
                      cmap = 'bone',
                      sdir = sdir, 
                      scale = ls[scale],
                      aspect = 'equal',
                     norm = None)
        
 
    @property
    def divergence(self):
        self.set_electric_field_representation('a')
        div_x, div_y = self.fwhm
        self.set_electric_field_representation('c')
        return (1.66*div_x)/2, (1.66*div_y)/2
    
    @property
    def pulse_energy(self):
        self.set_electric_field_representation('t')
        e = calc_pulse_energy(self)[0]
        self.set_electric_field_representation('f')
        return e
        
        
    def plot_phase(self, scale = "mm", label = "", title = "", context = 'talk', sdir = None):
        
        phase = self.get_phase().sum(-1)
        
        ls = {"m": 1,
              "cm": 1e2,
              "mm": 1e3,
              "um": 1e6,
              "nm": 1e9}
        
        mesh = get_mesh(phase, *self.get_spatial_resolution())

        colorbar_plot(phase, mesh,
                      label = label,
                      title = title,
                      xlabel = "x ({})".format(scale),
                      ylabel = "y ({})".format(scale),
                      clabel = "Phase",
                      vmin = -np.pi,
                      vmax = np.pi,
                      context = context,
                      cmap = 'hsv',
                      sdir = sdir, 
                      scale = ls[scale],
                      aspect = 'equal')
                      
    def get_beam_size(self, write = True, fraction = np.sqrt(np.log(1/2)/-2), threshold = 0.01, VERBOSE = False):
        
        px, py = self.get_spatial_resolution()
        ### fraction set to sqrt(ln(1/2)/-2) which equals the2sigma width for a guassian
        res, err = get_enclosed_energy(self.get_intensity().sum(-1), px, py, efraction = fraction, VERBOSE = VERBOSE,
                                       threshold = threshold)
       
        
        self.custom_fields['beam size'] = res
        if VERBOSE: print(self.custom_fields['beam size'])
        
        return res, err
 
    
    @property
    def x_axis(self):
        """
        return a one-dimensional array of wavefront coordinates along the horizontal-axis 
        """
        return np.linspace(self.xMin, self.xMax, self.nx)

    @property
    def y_axis(self):
        """
        return a one-dimensional array of wavefront coordinates along the vertical-axis 
        """
        return np.linspace(self.yMin, self.yMax, self.ny)

    @property
    def t_axis(self):
        """
        return a one-dimensional array of wavefront coordinates along the time-axis 
        """
        return np.linspace(self.zMin, self.zMax, self.nz)

    @property
    def f_axis(self):
        """
        return a one-dimensional array of wavefront coordinates along the
    temporal frequency-axis 
        """
        return 1/self.t_axis

 
    @property
    def get_mesh(self):
        return get_mesh(self.get_intensity().sum(-1), *self.get_spatial_resolution())
    
    def plot(self, scale = 'mm', sdir = None, label = None):
        
            ii = self.get_intensity().sum(-1)
            
            
            colorbar_plot(ii, mesh = self.get_mesh(),
                          label = label,
                  xlabel = "x ({})".format(scale),
                  ylabel = "y ({})".format(scale),
                  scale = ls[scale],
                  context = 'talk',
                  clabel = "Intensity (a.u.)",
                  grid = False,
                  aspect = 'equal',
                  vmax = np.max(ii),
                  sdir = sdir)

    def get_complex_radial_profile(self):
        """
        Calculate the radial profile of a complex array by azimuthal averaging:
            $$
            I_{radial}(R) = \int_0^R \frac{I(r)2\pi r}{\pi R^2} dr
            $$
        
        :param wfr: complex wavefield [np array]
        
        :returns prof: radial profile
        """
        wfr = self.as_complex_array()
        r = radial_profile(wfr[:,:,0].real, [wfr.shape[0]//2,wfr.shape[1]//2])[1]
        
        r = np.diag(r).copy()
        r[:r.shape[0]//2] *= -1
        
        rp = np.stack([radial_profile(wfr[:,:,i].real,
                                      [wfr.shape[0]//2,wfr.shape[1]//2])[0]
                       + radial_profile(wfr[:,:,i].imag,
                                        [wfr.shape[0]//2,wfr.shape[1]//2])[0]*1j
                       for i in range(wfr.shape[-1])])
        
        prof = np.moveaxis(rp, 0, -1)
        
        return prof, r
  


    def get_transverse_doc(self, VERBOSE = False):
        """
        get transverse degree of coherence of the wavefront across each of the
        transverse dimensions slices
        """
        
        p, r =  self.get_complex_radial_profile()
        nt = self.as_complex_array().shape[-1]
        J = np.dot(p, p.T.conjugate())/nt
        
        
        tdoc = np.diag(np.dot(J, J)).sum() / np.diag(J).sum()**2
        
        if VERBOSE:
            print("Transverse Degree of Coherence: {:.4f}".format(tdoc.real))
        
        self.custom_fields['tdoc'] = tdoc     
        return tdoc




    def get_coherence_len(self, VERBOSE = False):
        """
        Calculate coherence length of a complex wavefield of shape
        [nx, ny. nz]
        
        :param wfr: complex wavefield
        :param dx: horizontal pixel size
        :param dy: vertical pixel size
        
        :returns Jd: complex degree of coherence
        :returns clen: coherence length [m]
        """
        profile, r = self.get_complex_radial_profile()
        dx, dy = self.get_spatial_resolution()
        
        nt = self.as_complex_array().shape[-1]
        
        J = np.dot(profile, profile.T.conjugate())/ nt
        II = np.abs(np.diag(J))  # intensity as the main diagonal
        
        J /= II**0.5 * II[:, np.newaxis]**0.5
        Jd = np.abs(np.diag(np.fliplr(J)))  # DoC as the cross-diagonal
        
        lm = np.arange(Jd.shape[0])
    
        lm = lm[(lm >= Jd.shape[0]//2) & (Jd[lm] < 0.5)]
    
        rstep = np.sqrt((dx)**2 + (dy)**2)
    
        
        try:
            lm = lm[0] - Jd.shape[0]//2 
        except(IndexError):
            lm = np.inf
         
        clen = lm*rstep
    
        if VERBOSE: 
            print("Radial Coherence Length: {:.2f} um".format(clen*1e6))
        
        self.custom_fields['coherence length'] = clen
        
        return clen
    
    @property
    def extent(self):
        """
        returns matplotlib plt.imshow extent for plotting purposes.
        """
        return [self.params.Mesh.xMin, self.params.Mesh.xMax, self.params.Mesh.yMin, self.params.Mesh.yMax]
    
    @property
    def qxMin(self):
        return self.params.Mesh.qxMin
    
    @property
    def qyMin(self):
        return self.params.Mesh.qyMin
    
    @property
    def qxMax(self):
        return self.params.Mesh.qxMax
    
    @property
    def qyMax(self):
        return self.params.Mesh.qyMax
    
    @property
    def qx(self):
        return (self.qxMax-self.qxMin)/self.nx
    
    @property
    def qy(self):
        return (self.qyMax-self.qyMin)/self.ny
    
    
    @property
    def qx_axis(self):
        return np.linspace(self.qxMin, self.qxMax, self.nx)

    @property
    def qy_axis(self):
        return np.linspace(self.qyMin, self.qyMax, self.ny)
    
    @property
    def xMin(self):
        return self.params.Mesh.xMin
    
    @property
    def xMax(self):
        return self.params.Mesh.xMax
    
    @property
    def yMin(self):
        return self.params.Mesh.yMin
    
    @property
    def yMax(self):
        return self.params.Mesh.yMax    

    @property
    def dx(self):
        if self.params.wSpace == 'R-space':
            dx = (self.xMax-self.xMin)/self.nx
        elif self.params.wSpace == 'Q-space':
            dx = self.qx
        
        return dx
    
    @property
    def dy(self):
        if self.params.wSpace == 'R-space':
            dy = (self.yMax-self.yMin)/self.ny
        elif self.params.wSpace == 'Q-space':
            dy = self.qy
        
        return dy
    
    @property
    def nx(self):
        return self.params.Mesh.nx
    
    @property
    def ny(self):
        return self.params.Mesh.ny
    
    @property
    def nz(self):
        return self.params.Mesh.nSlices
   
    @property
    def x_axis(self):
        return np.linspace(self.extent[0], self.extent[1], self.nx)

    @property
    def y_axis(self):
        return np.linspace(self.extent[1], self.extent[2], self.ny) 
    
    @property
    def keys(self):
        return list(self.custom_fields.keys())
    
    @property
    def com(self):
        """ 
        return the center of mass, tuple - (comx, comy)
        """
        cx, cy = self.dx*self.nx//2, self.dy*self.ny//2
        
        return np.flip((get_com(self.get_intensity().sum(-1)
                        )-np.asarray([self.nx/2, self.ny//2])
                )*np.asarray([self.dx, self.dy]))
    
    @property
    def peak_intensity(self):
        return np.max(self.get_intensity().sum(-1))
    
    @property
    def source_properties(self):
        return self.custom_fields['source_properties']
    
    @property
    def metadata(self):
        return self.custom_fields['metadata']
    
    def get_y_profile(self, method = 'com'):
        """ 
        return a one-dimensional line profile along the x-axis
        
        :param method: determines how the line-profile is determined    
        """
        
        if method == 'com':
            x_profile = self.get_intensity().sum(-1)[:, int(self.nx//2 + self.com[0]//self.dx)]
        elif method == 'avg':
            x_profile = self.get_intensity().sum(-1).mean(1)
        elif method == 'c':
            x_profile = self.get_intensity().sum(-1)[:, self.nx//2]
        
        return x_profile
            
            
    
    def get_x_profile(self, method = 'com'):
        """ 
        return a one-dimensional line profile along the x-axis
        
        :param method: determines how the line-profile is determined    
        """
        
        if method == 'com':
            y_profile = self.get_intensity().sum(-1)[int(self.ny//2 + self.com[1]//self.dy),:]
        elif method == 'avg':
            y_profile = self.get_intensity().sum(-1).mean(0)
        elif method == 'c':
            y_profile = self.get_intensity().sum(-1)[self.ny//2, :]
        
        return y_profile
    
    @property        
    def fwhm(self):
        """
        return the fwhm of the intensity distribution by fitting a wavefront
        """
        
        if self.params.wSpace == 'R-space':
            x = self.x_axis
            y = self.y_axis
        
        if self.params.wSpace == 'Q-space':
            x = self.qx_axis
            y = self.qy_axis
        
        ix = self.get_x_profile()
        iy = self.get_y_profile()


        imax = self.peak_intensity       
        
        ### fit x
        initial_guess = [self.get_fwhm()[0], self.com[0], imax]
        p, co = fit_gaussian(x, ix, p0 = initial_guess, VERBOSE = False)
        
        fwhm_x = p[0]/sigma_to_fwhm
            
        ### fit y
        initial_guess = [self.get_fwhm()[1], self.com[1], imax]
        p, co = fit_gaussian(y, iy, p0 = initial_guess, VERBOSE = False)
        
        fwhm_y = p[0]/sigma_to_fwhm
        
        return abs(fwhm_x), abs(fwhm_y)    
    
    def gaussian_fit(self, method = 'com', VERBOSE = True):
        """
        return the properties of a guassian distribution fitted to the
        intensity distribution of the wavefront
        """
        
        x = self.x_axis
        y = self.y_axis
    
        ix = self.get_x_profile(method)
        iy = self.get_y_profile(method)
        
        imax = self.peak_intensity       
        
        ### fit x
        initial_guess = [self.get_fwhm()[0], self.com[0], imax]
        px, cox = fit_gaussian(x, ix, p0 = initial_guess, VERBOSE = VERBOSE)
        
        ### fit y
        initial_guess = [self.get_fwhm()[1], self.com[1], imax]
        py, coy = fit_gaussian(y, iy, p0 = initial_guess, VERBOSE = VERBOSE)


        return px, py
    

    def scale_beam_energy(self,  beam_energy = 1):
        self.data.arrEhor /= np.sqrt(self.pulse_energy/beam_energy)


if __name__ == '__main__':

# =============================================================================
#     from felpy.model.src.coherent import construct_SA1_pulse
#    
#     wfr = construct_SA1_pulse(200,200,4,1,.1)
#     wfr.metadata['trey'] = True
#     wfr.store_hdf5("../data/tmp/test.h5")
#     wfr.load_hdf5("../data/tmp/test.h5")
#     
# =============================================================================
    pass