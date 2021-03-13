# -*- coding: utf-8 -*-
import os
from wpg.beamline import Beamline as WPG_Beamline
from wpg import srwlib
from wpg import srwlpy as srwl

class Beamline(WPG_Beamline):
    
    def __init__(self, _srwl_wf = None):
        super().__init__(_srwl_wf)
        

    def propagate_sequential(self, wfr, outdir = None):
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
            
        for itr in range(len(self.propagation_options[0]['optical_elements'])):
            oe = self.propagation_options[0]['optical_elements'][itr]
            pp = self.propagation_options[0]['propagation_parameters'][itr]
            
            bl = srwlib.SRWLOptC([oe],[pp])
            
            try:
                print(oe.name)
            except:
                print(oe)
                
            srwl.PropagElecField(wfr._srwl_wf, bl)
            
 
            if outdir is None:
                wfr.view()


if __name__ == '__main__':

    from felpy.model.src.coherent import construct_SA1_pulse
    from felpy.model.tools import propagation_parameters
    from wpg.optical_elements import Drift
    wfr = construct_SA1_pulse(100,100,2,1,1)
    bl = Beamline()
    bl.append(Drift(10), propagation_parameters(1,1,1,1,'quadratic'))
    bl.propagate_sequential(wfr)