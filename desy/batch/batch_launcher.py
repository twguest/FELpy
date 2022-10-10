import inspect
import os 
import sys
import ast 
import json

from datetime import datetime

from felpy.model.beamlines.exfel_spb.methods import get_beamline_object
from felpy.model.beamline import Beamline
from felpy.utils.job_utils import batch_launcher

from labwork.about import dCache, logs

def propagate_NVE(wfr_directory, options):
    
    options = ast.literal_eval(options)
    
    save_directory, focus, analysis, filename, datetime, function, description = options
    analysis = json.loads(analysis.lower()) 
    
    print("info")
    print("wavefront directory: {}".format(wfr_directory))
    print("save directory: {}".format(save_directory))
    print("focus (i.e. beamline option): {}".format(focus))
    print("analysis: {}".format(analysis))
    print("datetime: {}".format(datetime))
    print("function: {}".format(function))
    print("description: {}".format(description))
    
  
    wfr = Wavefront()
    wfr.load_hdf5(wdir)
    sdir = sdir + "/{}/".format(focus)
    
    bl = get_beamline_object(ekev = 4.96, options = focus, crop = crop)
    
    wfr.custom_fields['focus'] = focus
    wfr.custom_fields['input directory'] = wfr_directory
    wfr.custom_fields['datetime'] = datetime
    wfr.custom_fields['function'] = function
    wfr.custom_fields['description'] = description
    wfr.custom_fields['bl'] = bl.__str__
            
    bl.propagate(wfr)
    
    if analysis: 
        wfr.analysis()
        
    wfr.store_hdf5(sdir)

    
if __name__ == '__main__':
    
    job_name = "4.98keV_Source_to_NVE"
    python_command = propagate_NVE
    input_directory = dCache + "/NanoKB-Pulse/source/"
    
    save_directory = dCache + "/NanoKB-Pulse/NVE/"
    mkdir_p(save_directory)
    
    log_directory = logs
    
    focus = "nano"
    analysis = True
    
    filename = __file__
    datetime = datetime.now().__str__()
    function = python_command.__name__
    
    description = "Propagate NanoKB Pulses 4.98 keV, 250 pC from Undulator-exit to the NVE"
    
    options = [save_directory,
              focus,
              analysis,
              filename,
              datetime,
              function,
              description]

    batch_launcher(python_command,
                   input_directory,
                   options = options,
                   nodes = 2, 
                   job_name = job_name,
                   partition = 'exfel',
                   log_directory = log_directory,
                   VERBOSE = True)