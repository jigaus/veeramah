#!/usr/bin/env python
# -*- coding: ASCII -*-

import string
import os
from subprocess import Popen,PIPE
from sys import argv
import gzip
import multiprocessing as mp
### imports modules that will be used

file_list=argv[1]#'file_list'
outdirectory=argv[2]#'/gpfs/scratch/kveeramah/stickleback/Batch_1_NovaSeq_process'
nb_reads=int(argv[3])#1000000
nbthreads=int(argv[4])#28
### defines arguments, number of reads and threads

fileout=open(outdirectory+'/'+file_list+'.FINISHED','w')

file=open(file_list,'r')
data=file.read()
data=string.split(data,'\n')
if data[-1]=='':
    del(data[-1])

runs=[]
samp_dirs={}
for g in range(len(data)):
    if data[g][-12:]=='.R1.fastq.gz':  #NYGC # looking for last twelve, that ending --> change to appro.
        k=string.split(data[g],'/')
        directory=string.join(k[:-1],'/')
        samp_dir=k[-2]
        run_id=string.split(k[-1],'.R1.fastq.gz')[0]
        if samp_dir not in samp_dirs:
            samp_dirs[samp_dir]=[]
        samp_dirs[samp_dir].append(run_id)
        batch_name=string.split(samp_dir,'_')[1]+'_'+str(samp_dirs[samp_dir].index(run_id)+1)
        runs.append([directory,samp_dir,run_id,batch_name])
    elif data[g][-8:]=='_1.fq.gz':   #BGI
        k=string.split(data[g],'/')
        directory=string.join(k[:-1],'/')
        samp_dir=k[-2]
        run_id=string.split(k[-1],'_1.fq.gz')[0]
        if samp_dir not in samp_dirs:
            samp_dirs[samp_dir]=[]
        samp_dirs[samp_dir].append(run_id)
        batch_name=samp_dir+'_'+str(samp_dirs[samp_dir].index(run_id)+1)
        runs.append([directory,samp_dir,run_id,batch_name])



samp_dirs_list=samp_dirs.keys()
for g in range(len(samp_dirs_list)):
    try:
        os.mkdir(outdirectory+'/'+samp_dirs_list[g])
    except:
        ok=1


nb_process=len(runs)

batches=[]
for g in range(0,nb_process,nbthreads):
    batches.append(runs[g:g+nbthreads])


cwd = os.getcwd()


def split_fastq(x,splits,outdirectory,nb_reads,output):
    os.chdir(cwd)

    directory=splits[x][0]
    samp_dir=splits[x][1]
    run_id=splits[x][2]
    batch_name=splits[x][3]

    Popen.wait(Popen('split_fastq.py '+run_id+' '+directory+' '+outdirectory+'/'+samp_dir+'/'+batch_name+' '+str(nb_reads),shell=True))

    output.put('finished '+batch_name)


for g in range(len(batches)):
    nbthreads2=len(batches[g])

    ###queue for parallelism output
    output = mp.Queue()

    # Setup a list of processes
    processes = [mp.Process(target=split_fastq, args=(x,batches[g],outdirectory,nb_reads,output)) for x in range(nbthreads2)]

    # Run processes
    for p in processes:
        p.start()

    # Exit the completed processes
    for p in processes:
        p.join()

    results = [output.get() for p in processes]
    print(results)

    for gg in range(len(results)):
        fileout.write(results[gg]+'\n')

fileout.close()

