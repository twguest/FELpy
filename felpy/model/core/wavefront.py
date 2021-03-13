from wpg import srwlpy 
from wpg.wavefront import Wavefront as WPG_Wavefront
from matplotlib import pyplot as plt
import imageio
import numpy as np
from wpg.wpg_uti_wf import calculate_fwhm
from felpy.analysis.optics.scalar.get_energy_statistics import get_energy_statistics

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
        return (h*c)/(self.params.photonEnergy)


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
        
        wDomain = self.params.wDomain
        self.set_electric_field_representation('a') 
        
        sig_x, sig_y = calculate_fwhm(self)['fwhm_x'], calculate_fwhm(self)['fwhm_y']
        
        self.set_electric_field_representation(wDomain) 

        return sig_x, sig_y

    def get_energy_statistics(self, integrate = False, VERBOSE = False, mpi = True):
        
        
        ii = self.get_intensity()
        
        if integrate:
            ii.sum(-1)
        
        dx, dy = self.get_spatial_resolution()
    
        dt = self.get_temporal_resolution()
        
        if integrate:
            dt *= self.params.Mesh.nSlices
        
        ekev = self.params.photonEnergy/1000
        
        return get_energy_statistics(ii, dx, dy, dt, ekev, mpi = mpi, VERBOSE = VERBOSE).run()        
        
        


if __name__ == '__main__':
    from felpy.model.src.coherent import construct_SA1_pulse
    wfr = construct_SA1_pulse(50,50,5,1,1)
    energy_statistics = wfr.get_energy_statistics(mpi = False)