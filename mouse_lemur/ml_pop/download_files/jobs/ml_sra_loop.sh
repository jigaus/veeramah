#!/bin/bash

#SBATCH --job-name=getml_sra
#SBATCH --ntasks-per-node=28
#SBATCH --nodes=2
#SBATCH --time=48:00:00
#SBATCH -p hbm-long-96core
#SBATCH --output=%j.getml_sra.out
#SBATCH --error=%j.getml_sra.error

cd $SLURM_SUBMIT_DIR

i=$(ml_accession_list.txt)
for i in $i
    prefetch $i --max-size u -O /gpfs/scratch/nguyen6/"$i"
    fasterq-dump $i 
done