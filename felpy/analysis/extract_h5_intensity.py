 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 19:11:57 2020

@author: twguest
"""

import sys
import os
import numpy as np

from utils.os_utils import mkdir_p
from utils.job_utils import JobScheduler
 
from felpy.model.core.wavefront import Wavefront

indir = "/opt/FELpy/felpy/data/test_pulses/"
outdir = "/opt/FELpy/felpy/data/"

tmp_dir = outdir + "/tmp/"
mkdir_p(tmp_dir)

intensity_dir = tmp_dir + "/integrated_intensity/"
mkdir_p(intensity_dir)

complex_dir = tmp_dir + "/complex_wavefield/"
mkdir_p(complex_dir)



def launch():
    """
    This part launches the jobs that run in main 
    """
    
    cwd = os.getcwd()
    script = os.path.basename(__file__)
    
    
    js = JobScheduler(cwd + "/" + script, logDir = "../../logs/",
                      jobName = "extractIntensity", partition = 'exfel', nodes = 4, jobType = 'array',
                      jobArray = indir)
    
    js.run(test = True)
    
    

def extract_intensity(fname):
    
    wfr = Wavefront()
    wfr.load_hdf5(indir + fname)
    
   
    np.save(intensity_dir + fname, wfr.get_intensity().sum(-1))
    np.save(complex_dir + fname, wfr.as_complex_array().sum(-1))
    
def compile_intensity_data():
    """
    note, it would be nice if this could read from squeue
    """

    
    np.save(outdir + "integrated_intensity" ,np.dstack(np.load(intensity_dir + ii) for ii in os.listdir(intensity_dir)))
    np.save(outdir + "complex_wavefield" ,np.dstack(np.load(intensity_dir + ii) for ii in os.listdir(intensity_dir)))
    

if __name__ == '__main__':
    
    mkdir_p(outdir)

    try:
        extract_intensity(sys.argv[1])
    except:
        for fname in os.listdir("/opt/FELpy/felpy/data/test_pulses/"):
            print(fname)
            extract_intensity(fname)

    compile_intensity_data()