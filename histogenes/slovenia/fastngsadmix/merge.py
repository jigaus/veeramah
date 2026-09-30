#!/usr/bin/env python
# -*- coding: ASCII -*-

import string
import os
from subprocess import Popen,PIPE
from sys import argv
import multiprocessing as mp
import random


nbthreads=22

chromosomes=['X']
random.shuffle(chromosomes)

nb_process=len(chromosomes)

batches=[]
for g in range(0,nb_process,nbthreads):
	batches.append(chromosomes[g:g+nbthreads])
	
#1240K_Gr37.chr11.bed
def process_aDNA(x,chrom_targs,output):
	chrom=chrom_targs[x]

	Popen.wait(Popen('merge_vcf_pyvcf_dnv_2.py 1240K_Gr37.chr'+chrom+'.bed ajd_chr'+chrom+' ',shell=True))
	Popen.wait(Popen('bgzip -f ajd_chr'+chrom+'.vcf',shell=True))
	Popen.wait(Popen('tabix -p vcf ajd_chr'+chrom+'.vcf.gz',shell=True))
	Popen.wait(Popen('make_homo_plink_edKV.py ajd_chr'+chrom+' ',shell=True))
	Popen.wait(Popen('vcftools --gzvcf ajd_chr'+chrom+'.vcf.gz --remove-indels --out ajd_chr'+chrom+'.PL --BEAGLE-PL --chr '+chrom+' ',shell=True))
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
