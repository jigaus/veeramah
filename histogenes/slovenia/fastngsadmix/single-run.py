#!/usr/bin/env python
# -*- coding: ASCII -*-

import string
import os
from subprocess import Popen,PIPE
from sys import argv
import multiprocessing as mp

indivlist=argv[1]
reference=argv[2]
nbthreads=int(argv[3])

individuals=[]
inds=open(indivlist,'r')
for line in inds:
	line=line.strip('\n').strip('\r').strip(' ')
	if len(line)>=1:
		individuals.append(line)

for z in range(len(individuals)):
	count=0
	for y in range(0,1):
		isExist = os.path.exists(str(y)+'/'+individuals[z]+'.qopt')
		if isExist:
			count+=1
	if count==50:
		individuals[z]=''
		
while('' in individuals):
	individuals.remove('')

nb_process=len(individuals)
nbthreads=len(individuals)
batches=[]
for g in range(0,nb_process,nbthreads):
	batches.append(individuals[g:g+nbthreads])


def process_aDNA(x,chrom_targs,output):
	chrom=chrom_targs[x]
	for j in range(0,1):
		isExist = os.path.exists(str(j)+'/'+chrom+'.qopt')
		if not isExist:
			# Popen.wait(Popen('fastNGSadmix -maxiter 10000 -likes /gpfs/scratch/dvyas/fastNGSadmix/unpruned/PLs/'+str(chrom)+'.PL -fname /gpfs/scratch/dvyas/fastNGSadmix/make_fname/refPanel_1240K_good.txt -Nname /gpfs/scratch/dvyas/fastNGSadmix/make_fname/nInd_1240K_good.txt -whichPops all -out '+str(j)+'/'+str(chrom)+' ',shell=True))
			Popen.wait(Popen('fastNGSadmix -maxiter 1000000 -likes /gpfs/scratch/dvyas/fastNGSadmix/unpruned/PLs/'+str(chrom)+'.PL -fname /gpfs/scratch/dvyas/fastNGSadmix/make_fname/refPanel_'+reference+'.txt -Nname /gpfs/scratch/dvyas/fastNGSadmix/make_fname/nInd_'+reference+'.txt -whichPops all -out '+str(j)+'/'+str(chrom)+' ',shell=True))
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




