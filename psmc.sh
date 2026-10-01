## PSMC nonsense
# tutorial: https://github.com/henriquevf/Tutorials/blob/main/PSMC%20tutorial.md

# start with bam + ref to make consensus sequence (fastq)
bcftools mpileup -Ou -f leon.fna leon.bam | bcftools call -c | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > <output.fq.gz>

# create psmc input file from fastq 
fq2psmcfa -q20 <input.fq.gz> > <output.psmcfa>


