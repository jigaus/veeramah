#!/bin/bash

cd /gpfs/scratch/nguyen6/calling/ajd
# Provide path to directory and specific EXT of files e.g. *.txt
search_dir=/gpfs/scratch/nguyen6/calling/ajd/*.bam
# Provide path to directory to save logs of Python
output_dir=/gpfs/scratch/nguyen6/calling
# Loop through all files in the directory
for file in $search_dir
do
    # Call the Python script and pass the current file as an argument
    echo "$file"
    python GenoCaller_indent_WC_dnv_SM_p2.py $file X /gpfs/scratch/nguyen6/calling/ajd/GRCh37.fa 5 0 1
done