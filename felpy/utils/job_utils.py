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
import inspect
import os

import numpy as np 

from felpy.utils.os_utils import mkdir_p
import shutil



def batch_launcher(python_command,
                   input_directory,
                   options = [],
                   nodes = 2, 
                   job_name = None,
                   partition = 'exfel',
                   log_directory = "",
                   VERBOSE = True):
    
    """
    A general batch launcher that calls on JobScheduler to submit batch jobs to slurm.
    
    The only functionality this adds from job_scheduler is that this should launch from w/in the file
    w/ process running as __main__ (likely clumsy)
    
    Note: the outer level function/file that calls this method will need to take data
    from sys.argv's: sys.argv[1] = wfr_directory, sys.argv[2:] = options[:]
    
    :param input_directory: input directory to loop over [str] - __file__ should load files
    :param save_directory
    :param options: script specific options
    
    
    """

    
    if job_name is None:
        
        if type(python_command) == type(random_string):
            job_name = python_command.__name__
    
    
    js = JobScheduler(python_command, log_directory = log_directory,
                      job_name = job_name, partition = 'exfel', nodes = nodes, job_type = 'array',
                      job_array = input_directory, options = options)
    
    js.run(test = False)
 

def random_string(length):
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
    
    def __init__(self, python_command, job_name, log_directory,
                 partition = 'exfel', nodes = 1,
                 job_type = 'single',
                 job_array = None,
                 n_spawns = 1,
                 VERBOSE = True,
                 runtime = "14-00:00:00",
                 email = "trey.guest@desy.de",
                 mailtype = "ALL",
                 options = None, 
                 rundir = None):
        """
        
        :param job_type: options = spawn, array, single 
        """
        
        if type(python_command) == str:
            self.command = "script"
        elif type(python_command) == type(random_string):
            self.command = "method"
        
        self.python_command = python_command
        self.job_name = job_name 
        self.partition = partition
        self.nodes = nodes
        self.job_type = job_type
        self.job_array = job_array
        self.n_spawns = n_spawns
        
        self.VERBOSE = VERBOSE
        
        
        if self.VERBOSE:
            print("\nInitialising Job Scheduler\n")
        if rundir:
            self.rundir = rundir
        else:
            self.rundir = os.getcwd()
            
        self.runTime = runtime
        self.log_directory = log_directory + job_name + "/"
        self.jobDir = log_directory + "jobs/" + job_name + "/"
        self.outDir = log_directory + "out/" + job_name + "/"
        self.errDir = log_directory + "error/" + job_name + "/"
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
    
        print("Python File: {}".format(self.python_command))
        print("Global Job Name: {}".format(self.job_name))
        print("# Nodes: {}".format(self.nodes))
        print("Partition: {}".format(self.partition))
        print("Job Type: {}".format(self.job_type))
        print("Run Dir: {}".format(self.rundir))
        print("Job Dir: {}".format(self.jobDir))
        print("Log Dir: {}".format(self.log_directory))
        print("Output Dir: {}".format(self.outDir))
        print("Error Dir: {}".format(self.errDir))
        print("Sending {} Diagnostics To: {}".format(self.mailtype, self.email))
        print("Script Options: {}".format(self.options))
        
        
        
        
    
    def test(self):
        """
        execute a single instance of the python script
        """
        
        jobFile = self.jobDir + "Test_{}.job".format(self.job_name)
        if self.VERBOSE:
            print("\nGenerating Python Script Test \n")
            print("Python Test Job: {}Test_{}".format(self.jobDir, jobFile))
            
        if self.job_type == 'spawn':
            array_item = np.random.randint(1e5)
        
        elif self.job_type == 'array':
            if type(self.job_array) == str:
                array_item = os.listdir(self.job_array)[np.random.randint(len(os.listdir(self.job_array)))]

            else:
                array_item = self.job_array[np.random.randint(len(self.job_array))]
        
        
        
        with open(jobFile, "w+") as fh:
        
            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --partition={} \n".format(self.partition))
    
            fh.writelines("#SBATCH --job-name=Test_{}.job\n".format(self.job_name))
            fh.writelines("#SBATCH --chdir {} \n".format(self.rundir))
            fh.writelines("#SBATCH --nodes={}\n".format(self.nodes))
            fh.writelines("#SBATCH --output={}Test_{}.out\n".format(self.outDir, self.job_name))
            fh.writelines("#SBATCH --error={}Test_{}.err\n".format(self.errDir, self.job_name))
            fh.writelines("#SBATCH --time={}\n".format(self.runTime))
            fh.writelines("#SBATCH --mail-type={}\n".format(self.mailtype))
            fh.writelines("#SBATCH --mail-user={}\n".format(self.email))
            
        
            if self.job_type == 'array' and array_item != None:
                
                if self.command == 'method':
                    method_dir = inspect.getfile(self.python_command)
                    run_directory, filename = method_dir.rsplit("/",)
                    fh.writelines("cd {}; python3 -c 'from {} import {};{}({})'".format(run_directory,filename.split(".py")[0],self.python_command.__name__, array_item))
                elif self.command == 'script':
                    fh.writelines("python3 {} {}".format(self.python_command, array_item))
            else:
                fh.writelines("python3 {}".format(self.python_command))
            
                        
            if self.options:
                for o in self.options:
                    fh.writelines(" {}".format(str(o)))
        
        
        fh.close()
        os.system("sbatch {}Test_{}.job".format(self.jobDir, self.job_name))
    
    
    def write_jobs(self, job_name, array_item = None):
        """
        wrapper for write_jobs data
        """

        
        jobFile = os.path.join(self.jobDir,"{}.job".format(job_name))
        
        with open(jobFile, "w+") as fh:
        
            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --partition={} \n".format(self.partition))
            fh.writelines("#SBATCH --chdir {} \n".format(self.rundir))
            fh.writelines("#SBATCH --nodes={}\n".format(self.nodes))
            fh.writelines("#SBATCH --output={}{}.out\n".format(self.outDir, job_name))
            fh.writelines("#SBATCH --error={}{}.err\n".format(self.errDir, job_name))
            fh.writelines("#SBATCH --time={}\n".format(self.runTime))
            fh.writelines("#SBATCH --mail-type={}\n".format(self.mailtype))
            fh.writelines("#SBATCH --mail-user={}\n".format(self.email))
            
            
            if self.job_type == 'array' and array_item != None:
                if self.command == 'method':
                    method_dir = inspect.getfile(self.python_command)
                    print(method_dir)
                    run_directory, filename = method_dir.rsplit("/",1)
                    fh.writelines("cd {}; python3 -c 'from {} import {};{}()' {}".format(run_directory,filename.split(".py")[0],self.python_command.__name__,
                                                                                           self.python_command.__name__,array_item))
                elif self.command == 'script':  
                    fh.writelines("python3 {} {}".format(self.python_command, array_item))
            elif self.job_type == 'spawn' and array_item != None:
                fh.writelines("python3 {} {}".format(self.python_command, array_item))
            else:
                fh.writelines("python3 {}".format(self.python_command))

                if self.options:
                    for o in self.options:
                        fh.writelines(" {}".format(str(o)))
        
        fh.close()
    
            
    def build_scripts(self):
        
        if self.VERBOSE: 
            print("\nGenerating {} Scripts".format(self.job_type))
            
        if self.job_type == 'spawn':
            
            for itr in range(self.n_spawns):
                
                seed = np.random.randint(1e05)
 
                jName = self.job_name + "_" + random_string(8)
                
                self.write_jobs(jName, array_item = seed)

                if self.VERBOSE:
                    print("Building Job File: {}.job".format(jName))              
            
        elif self.job_type == 'single':
            
            jName = self.job_name
            self.write_jobs(jName)
            
            
            if self.VERBOSE:
                print("Building Job File: {}.job".format(jName))              
        
            
        elif self.job_type == 'array':
            
            if type(self.job_array) == str:
                
                for array_item in os.listdir(self.job_array):
                    jName = self.job_name + "_" + array_item
                    self.write_jobs(jName, self.job_array+array_item)
                    
            if type(self.job_array) == list:

                for array_item in self.job_array:
                    
                    if type(array_item) == str:

                        jName = self.job_name  + "_" + array_item
                        self.write_jobs(jName, array_item)
                    else:

                        jName = self.job_name + array_item.__name__
                        self.write_jobs(jName, array_item.__name__)

                            
                if self.VERBOSE:
                    print("Building Job File: {}.job".format(jName))              
            
    def runScripts(self):
        
        print("\nSubmitting Jobs")
        for job in os.listdir(self.jobDir):
        
            os.system("sbatch {}".format(self.jobDir + job))

    
    def run(self, test = False):
        
        if test:
            self.test()
            
        self.build_scripts()
        self.runScripts()

        
def moduleTest():
    
    print("Testing Job Scheduler Module\n")
    
    print("Testing Job Scheduler Single Mode\n")
    js = JobScheduler(python_command = "../tests/testFile.py",
                      job_name = "SingleScheduleTest",
                      log_directory = "../logs/",
                      job_type = "single",
                      options = ["a", 2, 4, [1,2,3]])
    js.run(test = False)
  
    
    print("Testing Job Scheduler Spawn Mode\n")
    js = JobScheduler(python_command = "../tests/testFile.py",
                      job_name = "SpawnScheduleTest",
                      log_directory = "../logs/",
                      job_type = "spawn",
                      n_spawns = 20,
                      options = ["a", 2, 4, [1,2,3]])
    js.run(test = False)
     
        
    print("Testing Job Scheduler Array Mode\n")
    js = JobScheduler(python_command = "../tests/testFile.py",
                      job_name = "ArrayScheduleTest",
                      log_directory = "../logs/",
                      job_type = "array",
                      job_array = [1,2,3,4,5],
                      options = ["a", 2, 4, [1,2,3]])
    js.run(test = False)
  
  
    
if __name__ == "__main__":
    moduleTest()
