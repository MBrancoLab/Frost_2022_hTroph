#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=2G
#$ -t 1-16

module load samtools
module load bedtools

bam=$(sed -n "${SGE_TASK_ID}p" bam_list.txt)
sample=$(echo $bam | sed s/-dedup.bam//)

##collate read pairs
samtools collate -o $sample-collated.bam $bam

##convert to bed
bedtools bamtobed -bedpe -i $sample-collated.bam > $sample.bed

##clean, sort and convert to bedgraph
awk '$1==$4 && $6-$2 < 1000 {print $0}' $sample.bed > $sample.clean.bed
cut -f 1,2,6 $sample.clean.bed | sort -k1,1 -k2,2n -k3,3n > $sample.fragments.bed
bedtools genomecov -bg -i $sample.fragments.bed -g /data/Blizard-BrancoLab/Genomes/GRCh38/hg38.chrom.sizes > $sample.fragments.bedgraph
