from wpg import srwlpy 
from wpg.wavefront import Wavefront as WPG_Wavefront
from matplotlib import pyplot as plt
import imageio
import numpy as np
from wpg.wpg_uti_wf import calculate_fwhm
from felpy.analysis.optics.scalar.get_energy_statistics import get_energy_statistics
from felpy.utils.vis_utils import double_colorbar_plot, colorbar_plot
import numpy as np
from wpg.wpg_uti_wf import calculate_fwhm
import seaborn as sns
from felpy.utils.np_utils import get_mesh
from felpy.analysis.optics.scalar.enclosed_energy import get_enclosed_energy
from scipy.constants import h, c, e
from felpy.analysis.optics.scalar.centroid import get_com
from felpy.analysis.optics.complex.coherence import get_coherence_time
from felpy.model.core.fresnel_propagator import frensel_propagator
from felpy.utils.vis_utils import colorbar_plot
from felpy.utils.np_utils import get_mesh


ls = {"m": 1,
      "cm": 1e2,
      "mm": 1e3,
      "um": 1e6,
      "nm": 1e9}

class Wavefront(WPG_Wavefront):
    
    def __init__(self, _srwl_wf = None):
        super().__init__(_srwl_wf)
        
    def save_tif(self, outdir):
        imgReal = self.get_intensity()
        imageio.imwrite(outdir + ".tif", imgReal)


    def get_spatial_resolution(self):
        
        px = (self.params.Mesh.xMax - self.params.Mesh.xMin) / self.params.Mesh.nx
        py = (self.params.Mesh.yMax - self.params.Mesh.yMin) / self.params.Mesh.ny
        
        return px, py
    
    def get_temporal_resolution(self):
        
        wDomain = self.params.wDomain
        
        if wDomain == 'time':
            delta_tau = (self.params.Mesh.sliceMax - self.params.Mesh.sliceMin) / self.params.Mesh.nSlices
        else:
            self.set_electric_field_representation('t')
            delta_tau = (self.params.Mesh.sliceMax - self.params.Mesh.sliceMin) / self.params.Mesh.nSlices
            self.set_electric_field_representation(wDomain)

        return delta_tau

    def write(self, outdir):
        output = open(outdir + ".txt", "w")
        output.write(self.__str__())
        output.close()
        
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


    def get_wavelength(self):
        return (h*c)/(e*self.params.photonEnergy)


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
        
        self.set_electric_field_representation('f') 

        return sig_x, sig_y
    
    def get_fwhm(self):
        sig_x, sig_y = calculate_fwhm(self)['fwhm_x'], calculate_fwhm(self)['fwhm_y']
        return sig_x, sig_y


    def get_energy_statistics(self, integrate = False, VERBOSE = False, mpi = True, write = False):
        
        
        ii = self.get_intensity()
        
        if integrate:
            ii = ii.sum(-1)
            
        
        dx, dy = self.get_spatial_resolution()
    
        dt = self.get_temporal_resolution()
        
        if integrate:
            dt *= self.params.Mesh.nSlices
        
        ekev = self.params.photonEnergy/1000
        
        results = get_energy_statistics(ii, dx, dy, dt, ekev, mpi = mpi, VERBOSE = VERBOSE).run()    
        
        if write:
            self.custom_fields['pulse energy'] = results[0]
            self.custom_fields['nphotons'] = results[1]
        else:
            return results 
     
        
    def get_pulse_duration(self):
        
        self.set_electric_field_representation('t')
        t = self.params.Mesh.nSlices*self.get_temporal_resolution()
        self.set_electric_field_representation('f')
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
                      aspect = 'equal')
        
 
        
        
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
                      
    def get_beam_size(self, write = True):
        
        if write: 
            self.custom_fields['beam size'] = get_enclosed_energy(self.get_intensity().sum(-1), *self.get_spatial_resolution())
        else:
            return get_enclosed_energy(self.get_intensity().sum(-1), *self.get_spatial_resolution())
    
    def get_com(self, longitudinal = False, write = False):
        
        if longitudinal: 
            ii = self.get_intensity()
        else:    
            ii = self.get_intensity().sum(-1)
        
        idx = get_com(ii)
        
        px, py = self.get_spatial_resolution()
        
        if write:
            self.custom_fields['com'] = [px*idx[1], py*idx[0]]
        else:
            return [px*idx[1], py*idx[0]]
        
    def get_profile_1d(self):
        """
        return 1d profiles along the center of each transverse axis.
        """
        
        
        ii = self.get_intensity().sum(-1)
        
        idx = self.params.Mesh.nx//2
        idy = self.params.Mesh.ny//2
                
        
        ix = ii[:, int(idx)]
        iy = ii[int(idy), :]
        
        return ix, iy

    def get_coherence_time(self, mpi = False, write = False):
    
        self.set_electric_field_representation('t')
        time_step = self.get_temporal_resolution()
        
        if write:
            get_coherence_time(self.as_complex_array(), time_step, mpi = mpi)
        else:
            get_coherence_time(self.as_complex_array(), time_step, mpi = mpi)
            
 
    
    def propagate(self, z, upsample = 1):
        """
        returns propagated wfr
        """

        px, py = self.get_spatial_resolution()
        wf = frensel_propagator(self.as_complex_array().sum(-1), px, py,
                           self.get_wavelength(), z)
        return wf
    
    
    def get_mesh(self):
        mesh = get_mesh(self.get_intensity().sum(-1), *self.get_spatial_resolution())
        return mesh
    
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
    
if __name__ == '__main__':
    
    from felpy.model.src.coherent import construct_SA1_pulse
    wfr = construct_SA1_pulse(200,200,2,1,1)
    
    wfr.get_com(write = True)
    wfr.get_beam_size(write = True)
    wfr.get_energy_statistics(integrate = True, mpi = False, write = True)
    print(["{} {}".format(item, wfr.custom_fields[item]) for item in wfr.custom_fields])
    ef = wfr.propagate(1) 
    plt.imshow(abs(ef)**2)
    wfr.get_coherence_time()
    wfr.plot()