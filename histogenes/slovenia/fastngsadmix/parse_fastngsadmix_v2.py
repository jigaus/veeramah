#!/usr/bin/env python
# -*- coding: ASCII -*-

from sys import argv
import numpy as np
import random
import os
from subprocess import Popen,PIPE

folder=argv[1]
start=str(argv[2])
stop=str(argv[3])

outfile=open(folder+'.out','w')
command = "for i in {"+start+".."+stop+"} ; do grep ' after 0 runs' ${i}/*log >> "+folder+"_in.log; done"
command2 = "sed -i 's/ after 0 runs!//g' "+folder+"_in.log"
command3 = "sed -i 's/.log:best like /\t/g' "+folder+"_in.log"
command4 = "sed -i 's/\//\t/g' "+folder+"_in.log"

Popen.wait(Popen(command+' ',shell=True))
Popen.wait(Popen(command2+' ',shell=True))
Popen.wait(Popen(command3+' ',shell=True))
Popen.wait(Popen(command4+' ',shell=True))

outfile.write(command+'\n'+command2+'\n'+command3+'\n'+command4+'\n')

file=open(folder+"_in.log",'r')
#outfile=open(folder+'.out','w')
individual_dic={}

Popen.wait(Popen('mkdir -p '+folder+'/',shell=True))

for line in file:
	line=line.strip('\n')
	k=line.split('\t')
	run=int(k[0])
	ind=k[1]
	like=float(k[2])
	
	if ind not in individual_dic:
		individual_dic[ind]=[]
	if len(individual_dic[ind])<>run:
		print('list not sorted')
		print(line+'\n')
		quit()
	individual_dic[ind].append(like)
	
for i in individual_dic:
	outfile.write(i+'\t')
	outfile.write(str(np.max(individual_dic[i]))+'\t')
	
	finding = np.where(individual_dic[i]==np.max(individual_dic[i]))[0]
	if len(finding)>1:
		index=int(random.choice(finding))
		outfile.write(str(index)+' (Multiple)'+'\n')
	else:
		index=int(np.where(individual_dic[i]==np.max(individual_dic[i]))[0])
		outfile.write(str(index)+'\n')
	Popen.wait(Popen('cp '+str(index)+'/'+str(i)+'.* '+folder+'/ ',shell=True))
