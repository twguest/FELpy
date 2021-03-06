#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "1.0.1"
__maintainer__ = "Trey Guest"
__email__ = "twguest@students.latrobe.edu.au"
__status__ = "Developement"
"""

    
import string
import random
 
import os

import numpy as np 

from .os_utils import mkdir_p

import shutil
# =============================================================================
# 
# def timeUtils: 
#         stime = time.time() ## start time for logging
#     
#         os.system("sbatch {}Test-{}".format(self.jobDir, self.jobName))
#         
#         ftime = time.time() ## finish time for logging
#         
#         tdiff = ftime-stime ## time difference (time per run)
#         
# =============================================================================





def randomString(length):
    """
    writes a random string with length characters
    """

    
    letters = string.ascii_letters
    result_str = ''.join(random.choice(letters) for i in range(length))
    return result_str



class JobScheduler:
    """
    A Python class for scheduling Slurm jobs
    """
    
    def __init__(self, pycmd, jobName, logDir,
                 partition = 'exfel', nodes = 1,
                 jobType = 'single',
                 jobArray = None,
                 nSpawn = 1,
                 VERBOSE = True,
                 runtime = "14-00:00:00",
                 email = "trey.guest@desy.de",
                 mailtype = "ALL",
                 options = None, 
                 rundir = None):
        """
        
        :param jobType: options = spawn, array, single 
        """
        
        
        self.pycmd = pycmd
        self.jobName = jobName 
        self.partition = partition
        self.nodes = nodes
        self.jobType = jobType
        self.jobArray = jobArray
        self.nSpawn = nSpawn
        self.VERBOSE = VERBOSE
        
        
        if self.VERBOSE:
            print("\nInitialising Job Scheduler\n")
        if rundir:
            self.rundir = rundir
        else:
            self.rundir = os.getcwd()
        self.runTime = runtime
        self.logDir = logDir + jobName + "/"
        self.jobDir = logDir + "jobs/" + jobName + "/"
        self.outDir = logDir + "out/" + jobName + "/"
        self.errDir = logDir + "error/" + jobName + "/"
        self.email = email
        self.mailtype = mailtype
        self.options = options
        
        
        if os.path.exists(self.jobDir):
            shutil.rmtree(self.jobDir)
        mkdir_p(self.jobDir)
        
        if os.path.exists(self.outDir):
            shutil.rmtree(self.outDir)
        mkdir_p(self.outDir)
        
        if os.path.exists(self.errDir):
            shutil.rmtree(self.errDir)
        mkdir_p(self.errDir)
        
        if self.VERBOSE == True:
            self.__str__()
    
    def __str__(self):
    
        print("Python File: {}".format(self.pycmd))
        print("Global Job Name: {}".format(self.jobName))
        print("# Nodes: {}".format(self.nodes))
        print("Partition: {}".format(self.partition))
        print("Job Type: {}".format(self.jobType))
        print("Run Dir: {}".format(self.rundir))
        print("Job Dir: {}".format(self.jobDir))
        print("Log Dir: {}".format(self.logDir))
        print("Output Dir: {}".format(self.outDir))
        print("Error Dir: {}".format(self.errDir))
        print("Sending {} Diagnostics To: {}".format(self.mailtype, self.email))
        print("Script Options: {}".format(self.options))
        
        
        
        
    
    def test(self):
        """
        execute a single instance of the python script
        """
        
        jobFile = self.jobDir + "Test_{}.job".format(self.jobName)
        print(jobFile)
        if self.VERBOSE:
            print("\nGenerating Python Script Test \n")
            print("Python Test Job: {}Test_{}".format(self.jobDir, jobFile))
            
        if self.jobType == 'spawn':
            arrItem = np.random.randint(1e5)
        
        elif self.jobType == 'array':
            if type(self.jobArray) == str:
                arrItem = os.listdir(self.jobArray)[np.random.randint(len(os.listdir(self.jobArray)))]

            else:
                arrItem = self.jobArray[np.random.randint(len(self.jobArray))]
        
        
        
        with open(jobFile, "w+") as fh:
        
            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --partition={} \n".format(self.partition))
    
            fh.writelines("#SBATCH --job-name=Test_{}.job\n".format(self.jobName))
            fh.writelines("#SBATCH --chdir {} \n".format(self.rundir))
            fh.writelines("#SBATCH --nodes={}\n".format(self.nodes))
            fh.writelines("#SBATCH --output={}Test_{}.out\n".format(self.outDir, self.jobName))
            fh.writelines("#SBATCH --error={}Test_{}.err\n".format(self.errDir, self.jobName))
            fh.writelines("#SBATCH --time={}\n".format(self.runTime))
            fh.writelines("#SBATCH --mail-type={}\n".format(self.mailtype))
            fh.writelines("#SBATCH --mail-user={}\n".format(self.email))
            
        
            if self.jobType == 'array' and arrItem != None:
                fh.writelines("python {} {}".format(self.pycmd, arrItem))
            else:
                fh.writelines("python {}".format(self.pycmd))
            
                        
            if self.options:
                for o in self.options:
                    fh.writelines(" {}".format(str(o)))
        
        
        fh.close()
        os.system("sbatch {}Test_{}.job".format(self.jobDir, self.jobName))
    
    def jobScript(self, jobName, arrItem = None):
        """
        wrapper for jobscript data
        """
        jobFile = os.path.join(self.jobDir,"{}.job".format(jobName))
        
        with open(jobFile, "w+") as fh:
        
            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --partition={} \n".format(self.partition))
    
            fh.writelines("#SBATCH --job-name={}.job\n".format(jobName))
            fh.writelines("#SBATCH --chdir {} \n".format(self.rundir))
            fh.writelines("#SBATCH --nodes={}\n".format(self.nodes))
            fh.writelines("#SBATCH --output={}{}.out\n".format(self.outDir, jobName))
            fh.writelines("#SBATCH --error={}{}.err\n".format(self.errDir, jobName))
            fh.writelines("#SBATCH --time={}\n".format(self.runTime))
            fh.writelines("#SBATCH --mail-type={}\n".format(self.mailtype))
            fh.writelines("#SBATCH --mail-user={}\n".format(self.email))
            
        
            if self.jobType == 'array' and arrItem != None:
                fh.writelines("python {} {}".format(self.pycmd, arrItem))
            elif self.jobType == 'spawn' and arrItem != None:
                fh.writelines("python {} {}".format(self.pycmd, arrItem))
            else:
                fh.writelines("python {}".format(self.pycmd))
            
            
            
            if self.options:
                for o in self.options:
                    fh.writelines(" {}".format(str(o)))
        
        fh.close()
     
            
    def buildScripts(self):
        
        if self.VERBOSE: 
            print("\nGenerating {} Scripts".format(self.jobType))
            
        if self.jobType == 'spawn':
            
            for itr in range(self.nSpawn):
                
                seed = np.random.randint(1e05)
                 
                
                jName = self.jobName + "_" + randomString(8)
                
                self.jobScript(jName, arrItem = seed)

                if self.VERBOSE:
                    print("Building Job File: {}.job".format(jName))              
            
        elif self.jobType == 'single':
            
            jName = self.jobName
            self.jobScript(jName)
            
            
            if self.VERBOSE:
                print("Building Job File: {}.job".format(jName))              
        
            
        elif self.jobType == 'array':
            
            if type(self.jobArray) == str:
                
                for arrItem in os.listdir(self.jobArray):
                    print(arrItem)
                    jName = self.jobName + arrItem
                    self.jobScript(jName, arrItem)
                    
            if type(self.jobArray) == list:
                
                for arrItem in self.jobArray:
                    
                    if type(arrItem) == str:
                        print(arrItem)
                        jName = self.jobName + arrItem
                        self.jobScript(jName, arrItem)
                    else:
                        print(arrItem)
                        jName = self.jobName + arrItem.__name__
                        self.jobScript(jName, arrItem.__name__)
                
            else:
                
                for arrItem in self.jobArray:
                    jName = self.jobName + "_" + str(arrItem)
                    self.jobScript(jName, arrItem)
                            
                if self.VERBOSE:
                    print("Building Job File: {}.job".format(jName))              
            
    def runScripts(self):
        
        print("\nSubmitting Jobs")
        for job in os.listdir(self.jobDir):
        
            os.system("sbatch {}".format(self.jobDir + job))

    
    def run(self, test = False):
        
        if test:
            self.test()
            
        self.buildScripts()
        self.runScripts()

        
def moduleTest():
    
    print("Testing Job Scheduler Module\n")
    
    print("Testing Job Scheduler Single Mode\n")
    js = JobScheduler(pycmd = "../tests/testFile.py",
                      jobName = "SingleScheduleTest",
                      logDir = "../logs/",
                      jobType = "single",
                      options = ["a", 2, 4, [1,2,3]])
    js.run(test = False)
  
    
    print("Testing Job Scheduler Spawn Mode\n")
    js = JobScheduler(pycmd = "../tests/testFile.py",
                      jobName = "SpawnScheduleTest",
                      logDir = "../logs/",
                      jobType = "spawn",
                      nSpawn = 20,
                      options = ["a", 2, 4, [1,2,3]])
    js.run(test = False)
     
        
    print("Testing Job Scheduler Array Mode\n")
    js = JobScheduler(pycmd = "../tests/testFile.py",
                      jobName = "ArrayScheduleTest",
                      logDir = "../logs/",
                      jobType = "array",
                      jobArray = [1,2,3,4,5],
                      options = ["a", 2, 4, [1,2,3]])
    js.run(test = False)
  
  
    
if __name__ == "__main__":
    moduleTest()
