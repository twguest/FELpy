 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 19:11:57 2020

@author: twguest
"""

###############################################################################
import sys
sys.path.append("/opt/WPG/") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/WPG") # DESY MAXWELL PATH

sys.path.append("/opt/spb_model") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/spb_model") # DESY MAXWELL PATH
###############################################################################
###############################################################################
import os
import numpy as np

from utils.os_utils import mkdir_p
from utils.job_utils import JobScheduler
 
from wpg.wavefront import Wavefront

indir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/out"
outdir1 = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/data/ii/"
outdir2 = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/data/ii/cmplx/"

mkdir_p(outdir1)
mkdir_p(outdir2)

def launch():
    """
    This part launches the jobs that run in main 
    """
    
    cwd = os.getcwd()
    script = os.path.basename(__file__)
    
    
    js = JobScheduler(cwd + "/" + script, logDir = "../../logs/",
                      jobName = "extractIntensity", partition = 'exfel', nodes = 1, jobType = 'array',
                      jobArray = indir)
    
    js.run(test = True)
    
    
def getIntensity(fname):
    
    wfr = Wavefront()
    wfr.load_hdf5(indir + fname)
    
   
    np.save(outdir1 + fname, wfr.get_intensity())
    np.save(outdir2 + fname, wfr.toComplex()[0,:,:,:])
    
if __name__ == '__main__':
    
    
    getIntensity(sys.argv[1])
    
    
