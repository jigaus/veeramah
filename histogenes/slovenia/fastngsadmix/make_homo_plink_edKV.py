#!/usr/bin/env python
# -*- coding: ASCII -*-

import string
import gzip
from random import randint
from sys import argv
import numpy as np
from subprocess import Popen,PIPE

stem=argv[1]

file=gzip.open(stem+'.vcf.gz','r')
data=file.read()
data=string.split(data,'\n')
if data[-1]=='':
	del(data[-1])

fileout=open(stem+'.tped','w')


print 'starting making plink'
for g in range(len(data)):
	k=string.split(data[g],'\t')
	if data[g][0]=='#':
		if data[g][1]=='C':
			samps=k[9:]
			out=''
			for gg in range(len(samps)):
				out=out+samps[gg]+'\t'+samps[gg]+'\t0\t0\t2\t0\n'
			fileout2=open(stem+'.tfam','w')
			fileout2.write(out)
			fileout2.close()

	else:
		out=k[0]+' '+k[0]+':'+k[1]+' 0 '+k[1]

		ref=k[3]
		alt=k[4]
		
		if len(ref) > 1 or len(alt) > 1:
			continue
		key=string.split(k[8],':')
		GT_key=key.index('GT')
			
		genos=k[9:]
		for gg in range(len(genos)):
			GT=string.split(genos[gg],':')[GT_key]
			if GT=='./.':
				out=out+' 0 0'
			else:
				if GT=='0/0':
					out=out+' '+ref+' '+ref
				elif GT=='1/1':
					out=out+' '+alt+' '+alt
				elif GT=='0/1':
					AD_key=key.index('AD')
					AD=string.split(genos[gg],':')[AD_key]
					AD1,AD2=string.split(AD,',')
					AD1=float(AD1)
					AD2=float(AD2)
					if AD2 > AD1:
						out=out+' '+alt+' '+alt
					elif AD2 < AD1:
						out=out+' '+ref+' '+ref
					else:
						test=randint(0,1)
						if test==0:
							out=out+' '+ref+' '+ref
						elif test==1:
							out=out+' '+alt+' '+alt
				
				
		fileout.write(out+'\n')


fileout.close()


#Popen.wait(Popen('plink2 --tfile '+stem+' --make-bed --out '+stem,shell=True))

#Popen.wait(Popen('rm '+stem+'.tped',shell=True))
#Popen.wait(Popen('rm '+stem+'.tfam',shell=True))

#Popen.wait(Popen('mv '+stem+'.bed ./plink/',shell=True))
#Popen.wait(Popen('mv '+stem+'.bim ./plink/',shell=True))
#Popen.wait(Popen('mv '+stem+'.fam ./plink/',shell=True))
				
		
		
