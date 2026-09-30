#!/usr/bin/env python
# -*- coding: ASCII -*-

import string
import os
from subprocess import Popen,PIPE
from sys import argv
import multiprocessing as mp
import random
import gzip

beaglename=argv[1]#'All_from_WGS_rev_allchr.PL.BEAGLE.PL.gz'#
# basefile=argv[2]#'basefile'#
# stem=argv[3]#'test'#
# indiv=argv[2]

# needed=['LJ-Gosp_3011','LJ-Gosp_3018','LJ-Gosp_3021','LJ-Gosp_3022','LJ-Gosp_3028','LJ-Gosp_3029','LJ-Gosp_3050','LJ-Gosp_3062','LJ-Gosp_3071-A','LJ-Gosp_3071-B','LJ-Gosp_3110','LJ-Gosp_3118','LJ-Gosp_3125-A','LJ-Gosp_3125-C','LJ-Gosp_3125-D','LJ-Gosp_3157-A','LJ-Gosp_3157-C','LJ-Gosp_3157','LJ-Gosp_3174','LJ-Gosp_3184','LJ-Gosp_3186','LJ-Gosp_3187','LJ-Gosp_3190','LJ-Gosp_3192','LJ-Gosp_3193','LJ-Gosp_8004','LJ-Gosp_8032','LJ-Gosp_8045','LJ-Gosp_8073','LJ-Gosp_8083','LJ-Gosp_8105','LJ-Gosp_8111-A','LJ-Gosp_8111-B','LJ-Gosp_8111-C','LJ-Gosp_8113','LJ-Gosp_8121','LJ-Gosp_8128','LJ-Gosp_8139','LJ-Gosp_8143','LJ-Gosp_8144','LJ-Gosp_8146']
john=True
chrom_targs=[]
orderdict={}
basedict={}

input=gzip.open(beaglename,'rt')
header=input.readline()
input.close()
headerlist=header.strip('\n').split('\t')
iidlist=[]
for i in range(4,len(headerlist),3):
	iidlist.append(headerlist[i])
max=0

# baselist=[]
# basedata=open(basefile,'r')
# for a in basedata:
	# a=a.strip('\n')
	# baselist.append(a)
# basedata.close()

#header=open('header_u','r')
# base=open(stem,'w')
basecolumns='1-3,'
number=1
for b in range(len(iidlist)): #header:
	iid=iidlist[b] #b.strip('\n')
	orderdict[iid]=number
		# base.write(iid+'\n')
	start=3+(number)*3-2
	stop=start+2
	basedict[iid]=str(start)+'-'+str(stop)
	basecolumns='1-3,'+str(start)+'-'+str(stop)
	if john:
		command='zcat '+beaglename+' | cut -f'+basecolumns+' > '+iid+'.PL '
		print(command)
		chrom_targs.append(command)
	# Popen.wait(Popen(command+' ',shell=True))
	
	number+=1
nbthreads=50



nb_process=len(chrom_targs)

batches=[]
for g in range(0,nb_process,nbthreads):
	batches.append(chrom_targs[g:g+nbthreads])
	
#1240K_Gr37.chr11.bed
def process_aDNA(x,chrom_targs,output):
	chrom=chrom_targs[x]

	Popen.wait(Popen(chrom+' ',shell=True))
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
