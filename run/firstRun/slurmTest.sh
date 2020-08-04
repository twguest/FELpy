#!/bin/bash
#SBATCH --partition=exfel
#SBATCH --time=14-00:00:00                           # Maximum time requested
#SBATCH --nodes=16                                 # Number of nodes
##### note: from SLurm news file 17.11.0rc1:
##### Change --workdir in sbatch to be --chdir as in all other commands (salloc, srun)
#####SBATCH --workdir   /gpfs/exfel/data/user/guestt/spb_model/run/firstRun    # directory must already exist!
#SBATCH --chdir   /gpfs/exfel/data/user/guestt/spb_model/run/firstRun      # directory must already exist!
#SBATCH --job-name  guestt
#SBATCH --output    guestt-%N-%j.out            # File to which STDOUT will be written
#SBATCH --error     guestt-%N-%j.err            # File to which STDERR will be written
#SBATCH --mail-type END                           # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user twguest@students.latrobe.edu.au       # Email to which notifications will be sent. It defaults to <userid@mail.desy.de> if none is set.

python /gpfs/exfel/data/user/guestt/spb_model/run/firstRun/nanoKB-4.96keV-test.py
