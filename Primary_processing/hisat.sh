#!/bin/sh
#$ -cwd
#$ -pe smp 4
#$ -l h_rt=12:0:0
#$ -l h_vmem=2G
#$ -t 1-3

module load use.own
module load hisat2
module load trimgalore
module load samtools

R1=$(sed -n "${SGE_TASK_ID}p" rna_list.txt)
base=$(echo $R1 | sed s/.fastq.gz//)

##trim
trim_galore --dont_gzip $R1

##align
hisat2 --no-softclip --no-unal -x ~/genomes/GRCh38/hisat/genome -U ${base}_trimmed.fq -S $base.sam
#hisat2 --no-softclip --no-unal -x ~/genomes/rheMac10/hisat/rheMac10 -U ${base}_trimmed.fq -S $base.sam
#hisat2 --no-softclip --no-unal -x ~/genomes/mm10/hisat/genome -U ${base}_trimmed.fq -S $base.sam

##filter based on MAPQ
samtools view -b -q 2 $base.sam > $base.bam

