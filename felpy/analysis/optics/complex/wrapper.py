# -*- coding: utf-8 -*-
"""
wraps some felpy code to handle wpg input 
"""
from wpg.srwlib import srwl 
 
from felpy.analysis.optics.complex.coherence import get_coherence_time as gct

def get_coherence_time(wfr, VERBOSE = True):
    
    srwl.SetRepresElecField(wfr._srwl_wf, 't')
    time_step = (wfr.params.Mesh.sliceMax - wfr.params.Mesh.sliceMin)/wfr.params.Mesh.nSlices
    return gct(wfr.toComplex(), time_step)


    