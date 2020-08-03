#!/bin/bash
#SBATCH --partition=exfel
#SBATCH --time=3-12:00:00      # request 1 day and 12 hours
#SBATCH --mail-type=END,FAIL   # send mail when the job has finished or failed
#SBATCH --nodes=16              # number of nodes
#SBATCH --output=%x-%N-%j.out  # per default slurm writes output to slurm-<jobid>.out. There are a number of options to customize the job
cd /gpfs/exfel/data/user/guestt/spb_model/run/firstRun/ # the actual script.
python nanoKB-4.96keV.py
