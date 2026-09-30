#!/bin/bash

#SBATCH --job-name=split
#SBATCH --ntasks-per-node=28
#SBATCH --nodes=2
#SBATCH --time=48:00:00
#SBATCH -p long-40core
#SBATCH --output=%j.split.out
#SBATCH --error=%j.split.error

module purge
module load anaconda/2

cd $SLURM_SUBMIT_DIR

chmod -x split_fastq_MP.py

python ./split_fastq_MP.py ml_seq_list.txt /gpfs/scratch/nguyen6/splitout/ 1000000 56
