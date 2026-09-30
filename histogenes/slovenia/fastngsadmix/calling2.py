#!/usr/bin/env python
# -*- coding: ASCII -*-

import string
import os
from subprocess import Popen,PIPE
from sys import argv
import gzip
import multiprocessing as mp
import random
# samplist='samples'#argv[1] #samp_list
# ref_genome='/gpfs/projects/VeeramahGroup/ref_genomes/GRCh37/GRCh37.fa'
# read_len='76' #str(argv[3]) #51
# nbthreads=int(argv[1])

samp=argv[1] #samp_list
genome=argv[2]
udg=argv[3]
samps=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']


if genome=='Hg19':
	ref_genome='/gpfs/projects/VeeramahGroup/ref_genomes/hg19/hg19_Krause.fa'
else:
	ref_genome='/gpfs/projects/VeeramahGroup/ref_genomes/GRCh37/GRCh37.fa'
nb_process=len(samps)
nbthreads=len(samps)
batches=[]
for g in range(0,nb_process,nbthreads):
	batches.append(samps[g:g+nbthreads])


cwd = os.getcwd()

# sftp://dvyas@milan.seawulf.stonybrook.edu/gpfs/scratch/dvyas/LHP3/Tap_437.A0101.TF1_rmdup.sorted.q30.bam

def aDNA_pipeline_Krause_PE(x,splits,output):
	chrom=samp
	i=splits[x] #[1]
	# isExist = os.path.exists(chrom+'.'+str(i)+'.indent5.emit_all.vcf.gz') 
	# if not isExist:
	Popen.wait(Popen('mkdir -p '+chrom+' ',shell=True))
	
	if udg=='T':
		isExist = os.path.exists(chrom+'.'+str(i)+'.indent5.emit_all.vcf.gz.tbi')
		isExist2 = os.path.exists(chrom+'.'+str(i)+'.indent5.emit_all.vcf') 		
		if isExist or isExist2:
			pass
		else:
			if genome=='Hg19':
				orig_i=i
				i='chr'+i			
			Popen.wait(Popen('GenoCaller_indent_WC_dnv_SM.py '+chrom+'.bam '+str(i)+' '+ref_genome+' 5 0 1 '+chrom+' ',shell=True))
			
			if genome=='Hg19':
				sed1="sed -i 's/^chr//g' "+chrom+"."+str(i)+".indent5.emit_all.vcf"
				# sed2="sed -i 's/^chr//g' "+chrom+"."+str(i)+".indent5.vcf"
				# sed3="sed -i 's/^chr//g' "+chrom+"."+str(i)+".indent5.haploid.emit_all.vcf_like"
				Popen.wait(Popen(sed1+' ',shell=True))
				# Popen.wait(Popen(sed2+' ',shell=True))
				# Popen.wait(Popen(sed3+' ',shell=True))

			Popen.wait(Popen('bgzip -f '+chrom+'.'+str(i)+'.indent5.emit_all.vcf',shell=True))
			# Popen.wait(Popen('bgzip -f '+chrom+'.'+str(i)+'.indent5.vcf',shell=True))
			# Popen.wait(Popen('gzip '+chrom+'.'+str(i)+'.indent5.haploid.emit_all.vcf_like',shell=True))
			Popen.wait(Popen('tabix -p vcf '+chrom+'.'+str(i)+'.indent5.emit_all.vcf.gz',shell=True))
			# Popen.wait(Popen('tabix -p vcf '+chrom+'.'+str(i)+'.indent5.vcf.gz',shell=True))
				
			if genome=='Hg19':
				Popen.wait(Popen('mv '+chrom+'.'+str(i)+'.indent5.emit_all.vcf.gz '+chrom+'.'+str(orig_i)+'.indent5.emit_all.vcf.gz',shell=True))
				# Popen.wait(Popen('mv '+chrom+'.'+str(i)+'.indent5.vcf.gz '+chrom+'.'+str(orig_i)+'.indent5.vcf.gz',shell=True))
				Popen.wait(Popen('mv '+chrom+'.'+str(i)+'.indent5.emit_all.vcf.gz.tbi '+chrom+'.'+str(orig_i)+'.indent5.emit_all.vcf.gz.tbi',shell=True))
				# Popen.wait(Popen('mv '+chrom+'.'+str(i)+'.indent5.vcf.gz.tbi '+chrom+'.'+str(orig_i)+'.indent5.vcf.gz.tbi',shell=True))
				# Popen.wait(Popen('mv '+chrom+'.'+str(i)+'.indent5.haploid.emit_all.vcf_like.gz '+chrom+'.'+str(orig_i)+'.indent5.haploid.emit_all.vcf_like.gz',shell=True))
	else:
		isExist = os.path.exists(chrom+'.'+str(i)+'.aDNA.emit_all.vcf.gz.tbi')
		isExist2 = os.path.exists(chrom+'.'+str(i)+'.aDNA.emit_all.vcf') 
		if isExist or isExist2:
			pass
		else:
			if genome=='Hg19':
				orig_i=i
				i='chr'+i	
			Popen.wait(Popen('aDNA_GenoCaller_WC_dnv_SM.py '+chrom+'.bam '+str(i)+' '+ref_genome+'  '+chrom+'_5pCtoT_freq.txt '+chrom+'_3pGtoA_freq.txt 0 1 '+chrom+' ',shell=True))
			
			if genome=='Hg19':
				sed1="sed -i 's/^chr//g' "+chrom+"."+str(i)+".aDNA.emit_all.vcf"
				# sed2="sed -i 's/^chr//g' "+chrom+"."+str(i)+".aDNA.vcf"
				# sed3="sed -i 's/^chr//g' "+chrom+"."+str(i)+".aDNA.haploid.emit_all.vcf_like"
				Popen.wait(Popen(sed1+' ',shell=True))
				# Popen.wait(Popen(sed2+' ',shell=True))
				# Popen.wait(Popen(sed3+' ',shell=True))

			Popen.wait(Popen('bgzip -f '+chrom+'.'+str(i)+'.aDNA.emit_all.vcf',shell=True))
			# Popen.wait(Popen('bgzip -f '+chrom+'.'+str(i)+'.aDNA.vcf',shell=True))
			# Popen.wait(Popen('gzip '+chrom+'.'+str(i)+'.aDNA.haploid.emit_all.vcf_like',shell=True))
			Popen.wait(Popen('tabix -p vcf '+chrom+'.'+str(i)+'.aDNA.emit_all.vcf.gz',shell=True))
			# Popen.wait(Popen('tabix -p vcf '+chrom+'.'+str(i)+'.aDNA.vcf.gz',shell=True))
			
			if genome=='Hg19':
				Popen.wait(Popen('mv '+chrom+'.'+str(i)+'.aDNA.emit_all.vcf.gz '+chrom+'.'+str(orig_i)+'.aDNA.emit_all.vcf.gz',shell=True))
				# Popen.wait(Popen('mv '+chrom+'.'+str(i)+'.aDNA.vcf.gz '+chrom+'.'+str(orig_i)+'.aDNA.vcf.gz',shell=True))
				Popen.wait(Popen('mv '+chrom+'.'+str(i)+'.aDNA.emit_all.vcf.gz.tbi '+chrom+'.'+str(orig_i)+'.aDNA.emit_all.vcf.gz.tbi',shell=True))
				# Popen.wait(Popen('mv '+chrom+'.'+str(i)+'.aDNA.vcf.gz.tbi '+chrom+'.'+str(orig_i)+'.aDNA.vcf.gz.tbi',shell=True))
				# Popen.wait(Popen('mv '+chrom+'.'+str(i)+'.aDNA.haploid.emit_all.vcf_like.gz '+chrom+'.'+str(orig_i)+'.aDNA.haploid.emit_all.vcf_like.gz',shell=True))

		
	# else:
		# pass
	output.put('finished '+str(i)+' '+str(chrom))


for g in range(len(batches)):
	nbthreads2=len(batches[g])

	###queue for parallelism output
	output = mp.Queue()

	# Setup a list of processes
	processes = [mp.Process(target=aDNA_pipeline_Krause_PE, args=(x,batches[g],output)) for x in range(nbthreads2)]

	# Run processes
	for p in processes:
		p.start()
	 
	# Exit the completed processes
	for p in processes:
		p.join()

	results = [output.get() for p in processes]
	print(results)
