#!/usr/bin/env python
# -*- coding: ASCII -*-

import string
import os
from subprocess import Popen,PIPE
from sys import argv
import multiprocessing as mp

chrom_targs=[]

infile=open('/gpfs/scratch/nguyen6/ml_fastq2/mlfastq2.txt','r')
for line in infile:
	line=line.strip('\n')
	chrom_targs.append(line)
	
print(chrom_targs)
nbthreads=96
chrom_targs.reverse()
nb_process=len(chrom_targs)

batches=[]
for g in range(0,nb_process,nbthreads):
	batches.append(chrom_targs[g:g+nbthreads])


# command='rclone copy /gpfs/scratch/dvyas/'+chrom+' SeaWulf_backup:seawulf_backup/'+chrom+' '
# Popen.wait(Popen(command+' ',shell=True))
# command='echo '+chrom+' >> /gpfs/scratch/dvyas/rclonedone/donefile'
# Popen.wait(Popen(command+' ',shell=True))
	
#1240K_Gr37.chr11.bed
def process_aDNA(x,chrom_targs,output):
	chrom=chrom_targs[x]
	#command="cd /gpfs/scratch/dvyas/"+chrom+" ;  find . -name \"*\" -exec touch {} \;"
	#Popen.wait(Popen(command+' ',shell=True))
	# command="touch /gpfs/scratch/dvyas/MP.py ; touch /gpfs/scratch/dvyas/MP.sb"
	# Popen.wait(Popen(command+' ',shell=True))
	
	if not os.path.exists(os.path.basename(chrom)):
		command='wget '+chrom+' '
		Popen.wait(Popen(command+' ',shell=True))
	# command='echo '+chrom+' >> /gpfs/scratch/dvyas/rclonedone/donefile'
	# Popen.wait(Popen(command+' ',shell=True))
	
	output.put('finished '+chrom)
	
	
for g in range(len(batches)):
	nbthreads2=len(batches[g])
	
	###queue for parallelism output
	output = mp.Queue()

	# Setup a list of processes
	processes = [mp.Process(target=process_aDNA, args=(x,batches[g],output)) for x in range(nbthreads2)]

	# Run processes
	for p in processes:
		p.start()
	 
	# Exit the completed processes
	for p in processes:
		p.join()

	results = [output.get() for p in processes]

	print(results)
	
	

