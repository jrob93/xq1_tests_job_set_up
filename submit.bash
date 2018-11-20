#!/bin/bash
#
# This is a simple submission script for a single job
#
# run on a single node
#SBATCH --nodes=1
# Just one task per node - all cores allocated 
#SBATCH --ntasks-per-node=1
# run on xq1
#SBATCH -p xq1
# give the job a name
#SBATCH -J jamiedemo000

# this submission script is in the same dir as the exe so just run it without paths

srun ./rebound

