import inspect
import os 
import sys
import ast 
import json

from datetime import datetime

from felpy.model.beamlines.exfel_spb.methods import get_beamline_object
from felpy.model.core.beamline import Beamline
from felpy.model.core.wavefront import Wavefront
from felpy.utils.job_utils import batch_launcher
from felpy.utils.os_utils import mkdir_p
from felpy.model.tools import propagation_parameters
from wpg.optical_elements import Drift
from labwork.about import dCache, logs

def propagate_NVE():
            
    wfr_directory = sys.argv[1].replace("*", "/")
 
    job_name = "NKB_4980eV_250pC_NVE_to_EHC"
    python_command = propagate_NVE
    input_directory = dCache + "/NanoKB-Pulse/NVE/"
    
    save_directory = input_directory.replace("/NVE/", "/EHC/")
    mkdir_p(save_directory)
    
    log_directory = logs
    
    focus = "nano"
    analysis = False
    
    filename = __file__
    dt = datetime.now().__str__()
    function = python_command.__name__
    
    description = "Propagate NanoKB Pulses 4.98 keV, 250 pC from Undulator NVE mirror to the EHC Screen"
    
    append = None
    
    print("info")
    print("wavefront directory: {}".format(wfr_directory))
    print("save directory: {}".format(save_directory))
    print("focus (i.e. beamline option): {}".format(focus))
    print("analysis: {}".format(analysis))
    print("datetime: {}".format(dt))
    print("filename: {}".format(filename))
    print("function: {}".format(function))
    print("description: {}".format(description))
    
  
    wfr = Wavefront()
    wfr.load_hdf5(wfr_directory)
     
    bl = Beamline()
    bl.append(Drift(2.2+3.5), propagation_parameters(1/3,1,1/3,1,'quadratic'))
    bl.propagate(wfr)
    
    wfr.custom_fields['focus'] = focus
    wfr.custom_fields['job name'] = job_name
    wfr.custom_fields['input directory'] = wfr_directory
    wfr.custom_fields['datetime'] = dt
    wfr.custom_fields['function'] = function
    wfr.custom_fields['filename'] = filename
    wfr.custom_fields['description'] = description
    #wfr.custom_fields['bl'] = bl.__str__
                
    if analysis: 
        wfr.analysis()
        
    wfr.store_hdf5(wfr_directory.replace("/NVE/", "/EHC/"))

    
if __name__ == '__main__':
    
    job_name = "NKB_4980eV_250pC_NVE_to_EHC"
    
    python_command = propagate_NVE
    
    input_directory = dCache + "/NanoKB-Pulse/NVE/"
 
    
    log_directory = logs
    
    focus = "nano"
    analysis = True
 
    ### run python_command w/ all files in input_directory
    batch_launcher(python_command,
                   input_directory,
                   options = [],
                   nodes = 8, 
                   job_name = job_name,
                   partition = 'exfel',
                   log_directory = log_directory,
                   VERBOSE = True)