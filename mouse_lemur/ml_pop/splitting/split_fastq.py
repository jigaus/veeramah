#!/usr/bin/env python
# -*- coding: ASCII -*-

import gzip
from subprocess import Popen,PIPE
from sys import argv

stem=argv[1]#'Mul18-1_TAAGGCGA-CTCTCTAT_AC6R3HANXX_L006_001'
directory=argv[2]#'/vault/veeramah/dog/Project_VEE_10797_B01_NAN_Lane.2015-04-30/Sample_Mul18-1/fastq'
fileout_stem=argv[3]
nb_reads=int(argv[4])

try:
        file1=gzip.open(directory+'/'+stem+'.R1.fastq.gz','r')
        file2=gzip.open(directory+'/'+stem+'.R2.fastq.gz','r')
except:
        file1=gzip.open(directory+'/'+stem+'_1.fq.gz','r')
        file2=gzip.open(directory+'/'+stem+'_2.fq.gz','r')

#fileout=open('./trim/'+stem+'_trim.fastq','w')



count_read=0
count_file=1

head1=file1.readline()
head2=file2.readline()

print 'Writing split '+str(count_file)
fileout1=open(fileout_stem+'_split_'+str(count_file)+'_R1.fq','w')
fileout2=open(fileout_stem+'_split_'+str(count_file)+'_R2.fq','w')

while head1<>'':
    seq1=file1.readline()
    seq2=file2.readline()
    line1=file1.readline()
    line2=file2.readline()
    qual1=file1.readline()
    qual2=file2.readline()
    out1=head1+seq1+line1+qual1
    out2=head2+seq2+line2+qual2
#    out=head1+seq1[:50]+'\n'+line1+qual1[:50]+'\n'
    fileout1.write(out1)
    fileout2.write(out2)
#    fileout.write(out)
    count_read+=1

    head1=file1.readline()
    head2=file2.readline()

    if count_read==nb_reads:
        fileout1.close()
        fileout2.close()
        Popen.wait(Popen('gzip '+fileout_stem+'_split_'+str(count_file)+'_R1.fq',shell=True))
        Popen.wait(Popen('gzip '+fileout_stem+'_split_'+str(count_file)+'_R2.fq',shell=True))

        count_file+=1
        print 'Writing split '+str(count_file)
        fileout1=open(fileout_stem+'_split_'+str(count_file)+'_R1.fq','w')
        fileout2=open(fileout_stem+'_split_'+str(count_file)+'_R2.fq','w')

        count_read=0


fileout1.close()
fileout2.close()
Popen.wait(Popen('gzip '+fileout_stem+'_split_'+str(count_file)+'_R1.fq',shell=True))
Popen.wait(Popen('gzip '+fileout_stem+'_split_'+str(count_file)+'_R2.fq',shell=True))

#fileout.close()
