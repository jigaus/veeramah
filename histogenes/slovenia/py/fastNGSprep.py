#!/usr/bin/env python
# -*- coding: ASCII -*-

from sys import argv
import numpy as np
import os
from subprocess import Popen,PIPE
import gzip

beaglename=argv[1]#'Procrustes_Fonyod1_prunedEUR_allchr.PL.BEAGLE.PL.gz'#

input=gzip.open(beaglename,'rt')
header=input.readline()
input.close()
headerlist=header.strip('\n').split('\t')
iidlist=[]
for i in range(4,len(headerlist),3):
	iidlist.append(headerlist[i])
max=0

number=1
for b in range(len(iidlist)): #header:
	iid=iidlist[b] #b.strip('\n')
	start=3+(number)*3-2
	stop=start+2
	basecolumns='1-3,'+str(start)+'-'+str(stop)+''
	
	command='zcat '+beaglename+' | cut -f'+basecolumns+' > '+iid(i)+'.PL '
	Popen.wait(Popen(command+' ',shell=True))
	command2="sed -i 's/\:/_/g' "+iid(i)+".PL "
	Popen.wait(Popen(command2+' ',shell=True))
	number+=1




