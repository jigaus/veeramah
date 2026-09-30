 #!/usr/bin/env python
# -*- coding: ASCII -*-

###This program merges individual emit all vcf files created by the aDNA_GenoCaller and GenoCaller_indent scripts for biallelic SNPs
###pysam needs to be installed. Two input files are needed.
##The first input file needs to be a tab delimited bed file containing the list of SNPs to be included in the output vcf, one per line, of the form <chrom><start><end><all1><all2><snp_name>
##e.g. 2\t136608645\t136608646\tG\tA\t2:136608646\n
##One of two alleles should be the reference allele, though the order does not matter
##The second input file contains a list of individual emit_all vcf files to be merged 
##e.g.
##../2072_LCT.RG.LCT.bed.aDNA.emit_all.vcf	samp_2072
##../84001_LCT.RG.LCT.bed.aDNA.emit_all.vcf	samp_84001
##../84005_LCT.RG.LCT.bed.aDNA.emit_all.vcf	samp_84005
###to run type:
##merge_vcf.py <SNP_bedfile> <list_of_emit_alls> <reference genome>
##Output will be a merged vcf file.

import gzip
import string
import pysam
import numpy as np
import math
from sys import argv
import vcf


snp_bedfile=argv[1] 
samp_list_file=argv[2] 
ref_file=argv[3] #argv[3]'/vault/public/GRCh37/GRCh37.fa'

ref_seq=pysam.FastaFile(ref_file)

def phred2prob(x):
	return 10.0**(-x/10.0)

def prob2phred(x):
	return -10*math.log10(x)

PL_snp_dic={}

#AA,AC,AG,AT,CC,CG,CT,GG,GT,TT

PL_snp_dic['AC']=np.array([0,1,4])
PL_snp_dic['AG']=np.array([0,2,7])
PL_snp_dic['AT']=np.array([0,3,9])
PL_snp_dic['CA']=np.array([4,1,0])
PL_snp_dic['CG']=np.array([4,5,7])
PL_snp_dic['CT']=np.array([4,6,9])
PL_snp_dic['GA']=np.array([7,2,0])
PL_snp_dic['GC']=np.array([7,5,4])
PL_snp_dic['GT']=np.array([7,8,9])
PL_snp_dic['TA']=np.array([9,3,0])
PL_snp_dic['TC']=np.array([9,6,4])
PL_snp_dic['TG']=np.array([9,8,7])

geno_list=['0/0','0/1','1/0','1/1']

geno_dic={}
geno_dic['0/0']=0
geno_dic['0/1']=1
geno_dic['1/0']=1
geno_dic['1/1']=2
geno_dic[0]='0/0'
geno_dic[1]='0/1'
geno_dic[2]='1/1'


file=open(samp_list_file,'r')
data=file.read()
data=string.split(data,'\n')
if data[-1]=='':
	del(data[-1])

files_use=[]
samps=[]

for g in range(len(data)): #Edited by DNV
	# file_use,samp=string.split(data[g],'\t')
	# files_use.append(file_use)
	# samps.append(samp)
	files_use.append(data[g])

####open vcf files in pyvcf
vcf_reader={}
for g in range(len(files_use)): #Edited by DNV
	vcf_reader[g]=vcf.Reader(open(files_use[g], 'r'))
	if len(vcf_reader[g].samples) > 1:
		print 'Too many samples in VCF '+str(g)
		quit()
	else:
		samps.append(vcf_reader[g].samples[0])

nb_samps=len(samps)

fileout=open(samp_list_file+'.vcf','w')


###set up output headers
header='##fileformat=VCFv4.1\n'
for i in range(len(ref_seq.lengths)):
	header=header+'##contig=<ID='+ref_seq.references[i]+',length='+str(ref_seq.lengths[i])+'>\n'
header=header+'##reference=file:'+ref_file+'\n'

header=header+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'

for g in range(len(samps)):
	header=header+'\t'+samps[g]

header=header+'\n'

fileout.write(header)


# ####open vcf files in pyvcf
# vcf_reader={}
# for g in range(len(files_use)):
	# vcf_reader[g]=vcf.Reader(open(files_use[g], 'r'))


####iterate through individual SNPs
file=open(snp_bedfile,'r')

data=file.read()
data=string.split(data,'\n')
if data[-1]=='':
	del(data[-1])


snp_count=0

for g in range(len(data)):
	k=string.split(data[g],'\t')
	chrom=k[0]
	start=int(k[1])
	end=int(k[2])
	all1=k[3]
	all2=k[4]
	if (all1<>'0') and (all2<>'0'):
		ref_all=ref_seq.fetch(k[0],start,end)
		if all1==ref_all:
			ref=all1
			alt=all2
		elif all2==ref_all:
			ref=all2
			alt=all1

		GT_array=np.zeros(nb_samps,dtype='int32')
		DP_array=np.zeros(nb_samps,dtype='int32')
		AD_array=np.zeros((nb_samps,2),dtype='int32')
		GQ_array=np.zeros(nb_samps,dtype='int32')
		PL_array=np.zeros((nb_samps,3),dtype='int32')

		GT_array[:]=-9
		DP_array[:]=-9
		AD_array[:]=-9
		GQ_array[:]=-9
		PL_array[:]=-9

		miss=0
		
		for gg in range(len(files_use)):
			try:
				for record in vcf_reader[gg].fetch(chrom,start,end):
					GT=record.samples[0]['GT']
					DP=record.samples[0]['DP']
					GQ=record.samples[0]['GQ']
					PL=np.array(record.samples[0]['PL'])
					try:
						AD=record.samples[0]['AD']
					except:
						AD=[DP,0]	
					if GT in geno_list:
						GT_array[gg]=geno_dic[GT]
						DP_array[gg]=DP
						if geno_dic[GT]>=0:
							if ((len(record.ALT)==1) and (record.ALT[0]==alt)) or (record.ALT[0]==None):
								AD_array[gg]=AD
								PL10=PL[PL_snp_dic[ref+alt]]
								if np.min(PL10)==0:
									PL_array[gg]=PL10
									GQ_array[gg]=GQ
								else:
	  
									prob=np.zeros(3)
									for ggg in range(len(PL10)):
										prob[ggg]=phred2prob(float(PL10[ggg]))
									PL_new=-10*np.log10(prob/np.max(prob))
									PL_new2=np.nan_to_num(PL_new)
									PL_array[gg]=int(PL_new2[0]),int(PL_new2[1]),int(PL_new2[2])
									GQ_array[gg]=np.sort(PL_array[gg])[1]
							
							elif (len(record.ALT)==1) and (record.ALT[0]<>alt): ##DNV added
								GT_array[gg]=-9	##DNV added				
							
			except:
				GT_array[gg]=-9
				DP_array[gg]=-9
				AD_array[gg]=-9
				GQ_array[gg]=-9
				PL_array[gg]=-9
				miss+=1
		AC=np.sum(GT_array[np.where(GT_array<>-9)[0]])
		out=chrom+'\t'+str(end)+'\t'+k[5]+'\t'+ref+'\t'+alt+'\t100\tPASS\t'
		out=out+'AC='+str(np.sum(GT_array[np.where(GT_array<>-9)[0]]))+':'
		if len(np.where(GT_array<>-9)[0]>0):
			out=out+'AF='+str(round(AC/(len(np.where(GT_array<>-9)[0])*2.0),2))
		else:
			out=out+'AF=0'
		#out=out+'\tGT:DP:GQ:PL'
		out=out+'\tGT:DP:AD:GQ:PL'

		for gg in range(len(GT_array)):
			if GT_array[gg]<>-9:
				#out=out+'\t'+geno_dic[GT_array[gg]]+':'+str(DP_array[gg])+':'+str(GQ_array[gg])+':'+str(PL_array[gg][0])+','+str(PL_array[gg][1])+','+str(PL_array[gg][2])
				out=out+'\t'+geno_dic[GT_array[gg]]+':'+str(DP_array[gg])+':'+str(AD_array[gg][0])+','+str(AD_array[gg][1])+':'+str(GQ_array[gg])+':'+str(PL_array[gg][0])+','+str(PL_array[gg][1])+','+str(PL_array[gg][2])

			else:
				out=out+'\t./.:.:.:.:.'

		fileout.write(out+'\n')
		
		snp_count+=1

		if snp_count%1000==0:
			print 'At position '+k[-1]
		
	else:
		print 'Problem with '+k[-1]+'. Skipping'
	
fileout.close()


					  
