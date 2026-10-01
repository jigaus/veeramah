#!/bin/bash

#SBATCH --job-name=getml_sra
#SBATCH --ntasks-per-node=28
#SBATCH --nodes=2
#SBATCH --time=12:00:00
#SBATCH -p long-40core
#SBATCH --output=%j.getml_sra.out
#SBATCH --error=%j.getml_sra.error

cd $SLURM_SUBMIT_DIR

module purge
module load anaconda/3
conda activate mouse_lemur

export PATH=$PWD/sratoolkit.3.4.1-ubuntu64/bin:$PATH

sed -i 's/\r$//' ml_accession_list.txt

i=$(cat ml_accession_list.txt)
for i in $i
do
	prefetch $i --max-size u -O /gpfs/scratch/nguyen6/ml_reads/"$i"
	fasterq-dump $i -O /gpfs/scratch/nguyen6/ml_reads/"$i"
done
