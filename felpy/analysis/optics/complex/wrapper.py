# -*- coding: utf-8 -*-
"""
wraps some felpy code to handle wpg input 
"""
from wpg.srwlib import srwl 
 
from felpy.analysis.optics.complex.coherence import get_coherence_time as gct
from felpy.analysis.optics.complex.coherence import get_coherence_len as gcl

def get_coherence_time(wfr, VERBOSE = True):
    
    srwl.SetRepresElecField(wfr._srwl_wf, 't')
    time_step = (wfr.params.Mesh.sliceMax - wfr.params.Mesh.sliceMin)/wfr.params.Mesh.nSlices
    return gct(wfr.toComplex(), time_step)

def get_coherence_len(wfr, VERBOSE = True):
    srwl.SetRepresElecField(wfr._srwl_wf, 't')
    return gcl(wfr.toComplex(), wfr.pixelsize()[0], wfr.pixelsize()[1])
    
if __name__ == '__main__':
    
    from felpy.model.src.coherent import construct_SA1_pulse
    
    wfr = construct_SA1_pulse(1024, 1024, 1, 5.0, 0.25)
    cl = get_coherence_len(wfr)
    