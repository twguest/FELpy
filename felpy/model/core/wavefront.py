from wpg.srw import srwlpy 
from wpg.wavefront import Wavefront as WPG_Wavefront
from matplotlib import pyplot as plt
import imageio
import numpy as np
from wpg.wpg_uti_wf import calculate_fwhm
from felpy.utils.vis_utils import double_colorbar_plot, colorbar_plot
import numpy as np
from wpg.wpg_uti_wf import calculate_fwhm, calc_pulse_energy
import seaborn as sns
from felpy.utils.np_utils import get_mesh
from felpy.analysis.optics.scalar.enclosed_energy import get_enclosed_energy
from scipy.constants import h, c, e
from felpy.analysis.optics.scalar.centroid import get_com
from felpy.analysis.optics.complex.coherence import get_coherence_time
from felpy.model.core.fresnel_propagator import frensel_propagator
from felpy.utils.vis_utils import colorbar_plot
from felpy.utils.np_utils import get_mesh
from felpy.model.tools import radial_profile
from datetime import datetime
import os



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
        imgReal = self.get_intensity()
        imageio.imwrite(outdir + ".tif", imgReal)
 
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
           
    def get_temporal_resolution(self, VERBOSE = False):
        
        wDomain = self.params.wDomain
        
        if wDomain == 'time':
            delta_tau = (self.params.Mesh.sliceMax - self.params.Mesh.sliceMin) / self.params.Mesh.nSlices
        else:
            self.set_electric_field_representation('t')
            delta_tau = (self.params.Mesh.sliceMax - self.params.Mesh.sliceMin) / self.params.Mesh.nSlices
            self.set_electric_field_representation(wDomain)


        self.custom_fields['temporal resolution'] = delta_tau
        
        if VERBOSE: self.custom_fields['temporal resolution']
        
        return delta_tau


    def to_txt(self, outdir):
        
  
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
        
        self.set_electric_field_representation('f') 
        
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
    
    def get_com(self, longitudinal = False, VERBOSE = False):
        
        if longitudinal: 
            ii = self.get_intensity()
        else:    
            ii = self.get_intensity().sum(-1)
        
        idx = get_com(ii)
        
        px, py = self.get_spatial_resolution()
        
        self.custom_fields['com'] = [px*idx[1], py*idx[0]]
    
        if VERBOSE: print(self.custom_fields['com'])
        
        return [px*idx[1], py*idx[0]]
    
        
    def get_profile_1d(self, method = sum):
        """
        return 1d profiles along the center of each transverse axis.
        """
        
        
        ii = self.get_intensity().method(-1)
        
        idx = self.params.Mesh.nx//2
        idy = self.params.Mesh.ny//2
                
        
        ix = ii[:, int(idx)]
        iy = ii[int(idy), :]
        
        return ix, iy

    def get_coherence_time(self, mpi = False, VERBOSE = False):
    
        self.set_electric_field_representation('t')
        time_step = self.get_temporal_resolution()
        
        ctime = get_coherence_time(self.as_complex_array(), time_step, mpi = mpi)
            
        self.custom_fields['coherence time'] = ctime
    
        if VERBOSE: print(self.custom_fields['coherence time'])
        return ctime
    
    
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

    def get_complex_radial_profile(self):
        """
        Calculate the radial profile of a complex array by azimuthal averaging:
        
            I_{radial}(R) = \int_0^R \frac{I(r)2\pi r}{\pi R^2} dr
        
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

    def get_extent(self):
        """
        returns matplotlib plt.imshow extent for plotting purposes.
        """
        return [self.params.Mesh.xMin, self.params.Mesh.xMax, self.params.Mesh.yMin, self.params.Mesh.yMax]
        
    def analysis(self, VERBOSE = False, DEBUG = False):
        """
        run a full analysis of the non-plotting utils and write to h5 file        
        """
        
        print("Running Full Analysis")
        if DEBUG:
            print("")
            print("**********************************************")
            print("Warning: Coherence Measurements are Disabled")
            print("**********************************************")
            print("")
        
        if DEBUG:
            VERBOSE = True
            
        if VERBOSE: print("PULSE DURATION")
        self.get_pulse_duration(VERBOSE = DEBUG)
        
        if VERBOSE: print("SPATIAL RESOLUTION")
        self.get_spatial_resolution(VERBOSE = DEBUG)

        if VERBOSE: print("TEMPORAL RESOLUTION")
        self.get_temporal_resolution(VERBOSE = DEBUG)

        if VERBOSE: print("WAVELENGTH")
        self.get_wavelength(VERBOSE = DEBUG)

        if VERBOSE: print("ENERGY STATISTICS")
        self.get_energy_statistics(VERBOSE = DEBUG)

        if VERBOSE: print("PULSE DURATION")
        self.get_pulse_duration(VERBOSE = DEBUG)

        if VERBOSE: print("BEAM SIZE")
        self.get_beam_size(VERBOSE = DEBUG)

        if VERBOSE: print("CENTER OF MASS")
        self.get_com(VERBOSE = DEBUG)

        if VERBOSE: print("COHERENCE TIME")
        ###self.get_coherence_time(VERBOSE = DEBUG, mpi = True)

        if VERBOSE: print("COHERENCE LENGTH")
        ###self.get_coherence_len(VERBOSE = DEBUG)

        if VERBOSE: print("TDOC")
        ###self.get_transverse_doc(VERBOSE = DEBUG)


        if DEBUG:
            print(self.custom_fields)
            
    
    def print(self):
        print(self.custom_fields)
        
    def get_keys(self):
        return list(self.custom_fields.keys())
    
    def get_values(self):
        return list(self.custom_fields.values())
    
    def log(self, bl = None, descriptor = None):
        
        self.custom_fields['user'] = os.getlogin()
        self.custom_fields['run directory'] = os.getcwd()
        self.custom_fields['datetime'] = datetime.now().__str__().split(".")[0]
        
        if bl is not None:
            self.custom_fields['beamline'] = bl.__str__()
            
        if descriptor is not None:
            self.custom_fields['function'] = descriptor[0]
            self.custom_fields['filename'] = descriptor[1] 
            self.custom_fields['summary'] = descriptor[2]
        
if __name__ == '__main__':
    
    pass 
    # from felpy.model.source.coherent import construct_SA1_pulse
   
    # wfr = construct_SA1_pulse(200,200,4,1,.1)
    # wfr.analysis()
 