#!/bin/sh
#$ -cwd
#$ -pe smp 4
#$ -l h_rt=4:0:0
#$ -l h_vmem=2G
#$ -t 1-4

module load use.own
module load hisat2
module load samtools
module load R
module load bedtools


R1=$(sed -n "${SGE_TASK_ID}p" rna_list.txt)
R2=$(echo $R1 | sed s/_R1/_R2/)
base=$(basename $R1 | sed s/_R1.fastq.gz//)


##align
hisat2 -x ~/genomes/GRCh38/hisat/genome -1 $R1 -2 $R2 -S $base.sam

##filter based on MAPQ
samtools view -b -q 2 $base.sam > $base.bam

##sort
samtools sort -@ 4 -o ${base}-sorted.bam $base.bam

##stringtie
~/stringtie/stringtie ${base}-sorted.bam --rf -G ~/genomes/GRCh38/gencode.v38.annotation.gtf -o $base.gtf

##get TSSs
Rscript get_tss.R $base.gtf ${base}-tss.bed

##intersect with Repeatmasker
intersectBed -wa -wb -a ~/genomes/GRCh38/Repeatmasker_hg38_4.0.5_noSimple.bed -b ${base}-tss.bed > ${base}-TEtranscripts.txt
