#!/bin/bash

#SBATCH --job-name=ajd
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH -p short-40core

# load modules
module purge
module load py/2.7.15

# Change to working directory
cd $SLURM_SUBMIT_DIR

# run
cd /gpfs/home/nguyen6/data
./merge.sh