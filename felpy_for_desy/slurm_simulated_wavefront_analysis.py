# -*- coding: utf-8 -*-
import sys
import os
from felpy.utils.job_utils import JobScheduler
from felpy.model.core.wavefront import Wavefront

from labwork.about import logs

PYCMD = os.getcwd() + "simulated_wavefront_analysis.py"



def launch_batch(in_directory, out_directory,
                 job_name,
                 logs = logs, 
                 nodes = 2):
    
    
    JobScheduler(pycmd = PYCMD,
                 jobName = job_name,
                 logDir = logs,
                 nodes = nodes,
                 jobType = 'array',
                 jobArray = os.listdir(in_directory),
                 options = out_directory)
    
def launch_single(in_directory, out_directory,
                 job_name,
                 logs = logs, 
                 nodes = 2):
   
    
    
    JobScheduler(pycmd = PYCMD,
                 jobName = job_name,
                 logDir = logs,
                 nodes = nodes,
                 jobType = 'array',
                 jobArray = in_directory,
                 options = out_directory)
    
   
if __name__ == '__main__':
    
    DEBUG = True
    
    if len(sys.argv) >= 2:
        in_directory = sys.argv[1]
    else:
        in_directory = input("Input Directory:")
        
    if len(sys.argv) >= 3:
        out_directory = sys.argv[2]
    else:
        out_directory = input("Output Directory:")
    
    if len(sys.argv) >= 4:
        job_name = sys.argv[3]
    else: 
        job_name = 'sim_wfr_analysis'
    
    
    if len(sys.argv) >= 5:
        logs = sys.argv[4]
    else: 
        logs = logs
    
    
    if len(sys.argv) >= 6:
        nodes = sys.argv[5]
    else: 
        nodes = 2
    
    if ".h5" or ".hdf5" in in_directory:
        mode = 'single'    
    else:
        mode = 'batch'

    if mode == 'single':
        
        launch_single(in_directory,
                      out_directory,
                      job_name,
                      logs,
                      nodes)
       
    elif mode == 'batch':
        launch_batch(in_directory,
              out_directory,
              job_name,
              logs,
              nodes)

    if DEBUG: 
        print(in_directory)
        print(out_directory)
        print(job_name)
        print(logs)
        print(nodes)