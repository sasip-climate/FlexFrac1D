#!/bin/bash

#OAR -l /cpu=1/core=2,walltime=24:00:00

# Modules loading

source /soft/env.bash

module load python/python3.9

# Launch compute job

python3 MF1DSpecExp.py > MF1D.out

# > /nfs_scratch/$USER/output_MF1D
