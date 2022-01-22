# -*- coding: utf-8 -*-
import sys
import os
from felpy.utils.job_utils import JobScheduler
from felpy.utils.os_utils import mkdir_p

from labwork.about import logs

PYCMD = os.getcwd() + "/propagate_EHC_data.py"



def launch_batch(in_directory, out_directory,
                 job_name,
                 logs = logs, 
                 nodes = 4):
    
    
    js = JobScheduler(pycmd = PYCMD,
                 jobName = job_name,
                 logDir = logs,
                 nodes = nodes,
                 jobType = 'array',
                 jobArray = os.listdir(in_directory),
                 options = out_directory)
    
    js.run(test = False)
    
def launch_single(in_directory, out_directory,
                 job_name,
                 logs = logs, 
                 nodes = 4):
   
    
    
    js = JobScheduler(pycmd = PYCMD,
                 jobName = job_name,
                 logDir = logs,
                 nodes = nodes,
                 jobType = 'single',
                 jobArray = in_directory,
                 options = [in_directory, out_directory])
    
    js.run(test = False)

   
if __name__ == '__main__':
    
    DEBUG = True
    
    if len(sys.argv) >= 2:
        in_directory = sys.argv[1]
    else:
        in_directory = input("Input Directory:")
        
    if len(sys.argv) >= 3:
        out_directory = sys.argv[2]
    else:
        out_directory = in_directory.replace("EHC", "PAM")
    
    if len(sys.argv) >= 4:
        job_name = sys.argv[3]
    else: 
        job_name = 'propagate_EHC_data'
    
    
    if len(sys.argv) >= 5:
        logs = sys.argv[4]
    else: 
        logs = logs
    
    
    if len(sys.argv) >= 6:
        nodes = sys.argv[5]
    else: 
        nodes = 2
    
    mode = input("batch of single?")
        
    if mode == 'single':
        print("Launching Single Job")
        launch_single(in_directory,
                      out_directory,
                      job_name,
                      logs,
                      nodes)
       
    elif mode == 'batch':
        print("Launching Batch Job")
        mkdir_p(out_directory)
        for item in os.listdir(in_directory):
            
            launch_single(in_directory + "{}".format(item),
                          out_directory + "{}".format(item),
                          job_name,
                          logs,
                          nodes)

    if DEBUG: 
        print(in_directory)
        print(out_directory)
        print(job_name)
        print(logs)
        print("nodes: ",nodes)