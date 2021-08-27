# -*- coding: utf-8 -*-
import os
from wpg.beamline import Beamline as WPG_Beamline
from wpg import srwlib
from wpg.srw import srwlpy as srwl
from wpg.wpg_uti_wf import plot_intensity_map

class Beamline(WPG_Beamline):
    
    def __init__(self, _srwl_wf = None):
        super().__init__(_srwl_wf)
        

    def propagate_sequential(self, wfr, outdir = None, return_intensity = False, return_mesh = False):
        """
        Propagate sequentially through each optical element in beamline.

        :param wfr: Input wavefront (will be re-writed after propagation)
        :param outdir: save directory
        """
        
        wfr.view()
        
        print(outdir)
        if outdir is not None:
            
            if os.path.exists(outdir) == False: 
                assert(ValueError, "Output Directory Could Not Be Found")
            else: 
                pass
        

        if outdir is not None:
            wfr.write(outdir + "initialSource")
            wfr.view()
        
        if return_intensity:
            intensity = []
        if return_mesh:
            mesh = []
            
        for itr in range(len(self.propagation_options[0]['optical_elements'])):
            oe = self.propagation_options[0]['optical_elements'][itr]
            pp = self.propagation_options[0]['propagation_parameters'][itr]
            
            bl = srwlib.SRWLOptC([oe],[pp])
            
            try:
                print(oe.name)
            except:
                print(oe)
                
            srwl.PropagElecField(wfr._srwl_wf, bl)
            
            if return_intensity:
                intensity.append(wfr.get_intensity().sum(-1))
            if return_mesh:
                mesh.append(wfr.get_mesh())
                    
 
            if outdir is None:
                if return_intensity:
                    pass
                else:
                    plot_intensity_map(wfr)

        if return_intensity:
            
            if return_mesh:
                return intensity, mesh
            else:
                return intensity
            
            
if __name__ == '__main__':

    from felpy.model.src.coherent import construct_SA1_pulse
    from felpy.model.tools import propagation_parameters
    from wpg.optical_elements import Drift
    wfr = construct_SA1_pulse(512,512,2,5.0,1)
    print(wfr.params.wDomain)
    wfr.plot()
    bl = Beamline()
    bl.append(Drift(10), propagation_parameters(2,1,2,1,'quadratic'))
    bl.propagate_sequential(wfr)