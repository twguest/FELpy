#!/usr/bin/env python

"""
Job Batch File Submission Script for NKB Zernike Fit
"""

import os

from os import listdir
indir = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/in/"

f = listdir(indir)

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)
    

job_directory = "/gpfs/exfel/data/user/guestt/spb_model/run/NanoKB-FastPulse/"
log_dir = "/gpfs/exfel/data/user/guestt/spb_model/logs/zernike/"

# Make top level directories
mkdir_p(job_directory)
mkdir_p(log_dir)

for fname in f:
    job_file = os.path.join(job_directory,"{}.job".format(fname))
    print("Launching: {} via {}".format(fname, job_file))

    with open(job_file, "w+") as fh:
        
        fh.writelines("#!/bin/bash\n")
        fh.writelines("#SBATCH --partition=exfel \n")

        fh.writelines("#SBATCH --job-name={}.job\n".format(fname))
        fh.writelines("#SBATCH --chdir {} \n".format(job_directory))
        fh.writelines("#SBATCH --nodes=1\n ")
        fh.writelines("#SBATCH --output={}{}.out\n".format(log_dir, fname))
        fh.writelines("#SBATCH --error={}{}.err\n".format(log_dir, fname))
        fh.writelines("#SBATCH --time=14-00:00:00\n")
        fh.writelines("#SBATCH --mail-type=ALL\n")
        fh.writelines("#SBATCH --mail-user=trey.guest@desy.de\n")
        fh.writelines("python /gpfs/exfel/data/user/guestt/spb_model/model/analysis/fitZernicke.py {}\n".format(fname) )
        
    fh.close()
    
    os.system("sbatch %s" %job_file)
