import inspect
import os 
import sys
import ast 
import json

from felpy.model.beamlines.exfel_spb.methods import get_beamline_object
from felpy.model.beamline import Beamline
from felpy.utils.job_utils import JobScheduler
from felpy.utils.os_utils import mkdir_p
from felpy.model.tools import propagation_parameters, scale
from felpy.model.wavefront import Wavefront

from wpg.wpg_uti_wf import plot_intensity_map as plt
from wpg.optical_elements import Drift, Aperture

from labwork.about import logs, dCache ### FOR PRIVATE USE, BOTH ARE DIRECTORIES


def propagate(wfr_directory, sdir, focus, analysis = False, crop = [], append = None, descriptor = "", VERBOSE = True):
    
    print("info")
    print("wavefront directory: {}".format(wfr_directory))
    print("save directory: {}".format(sdir))
    print("focus (i.e. beamline option): {}".format(focus))
    print("analysis: {}".format(analysis))
    print("crop: {}".format(crop))
    print("append: {}".format(append))
    
    for item in descriptor:
        print(item)
    
    wfr = Wavefront()
    wfr.load_hdf5(wdir)
    sdir = sdir + "/{}/".format(focus)
    
    bl = get_beamline_object(ekev = 4.96, options = focus, crop = crop)

    
    if append is not None:

        for item in append:
            bl.append(item[0],item[1])

    wfr.log(bl, descriptor)
    
    bl.propagate(wfr)
    
    if analysis: 
        wfr.analysis()
        
    wfr.store_hdf5(sdir)

    
    
def propagation_batch_launcher(input_directory, sdir, focus, analysis = False, crop = None, append = None, descriptor = "", VERBOSE = True):
    """
    This part launches the jobs that run in main 
    """
    
    wfr_directory = input_directory

    mkdir_p(sdir)
    
    cwd = os.getcwd()
    script =  __file__    
    filename = script
    function = inspect.currentframe().f_code.co_name
    
    
    js = JobScheduler(cwd + "/" + script, logDir = logs + "/",
                      jobName = "NVE_4.96keV", partition = 'exfel', nodes = 2, jobType = 'array',
                      jobArray = wfr_directory, options = [sdir, focus, analysis, crop, append, descriptor])
    
    js.run(test = False)
 
    
    
def FAST_NKB_NVE()
    
if __name__ == '__main__':
    
    wdir = sys.argv[1]
    sdir = sys.argv[2] 
    focus = sys.argv[3]
    
    analysis = json.loads(sys.argv[4].lower()) 
    
    print(sys.argv[5])
    if sys.argv[5] == "None":
        crop = None
    else:
        crop = ast.literal_eval(sys.argv[5])
    
    if sys.argv[6] == "None":
        append = None
    else:
        append = ast.literal_eval(sys.argv[6])
        
    description = ast.literal_eval(sys.argv[7])
    
    
    
    propagate(wdir, sdir, focus, analysis, crop, append, description)