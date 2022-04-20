#!/bin/sh
#$ -cwd
#$ -pe smp 4
#$ -l h_rt=48:0:0
#$ -l h_vmem=2G
#$ -t 1-29

module load bowtie2
module load trimgalore
module load samtools
module load python/3.6.3
module load use.own
module load localtools


R1=$(sed -n "${SGE_TASK_ID}p" fastq_list.txt)
R2=$(echo $R1 | sed s/_R1.fq.gz/_R2.fq.gz/)
base_name=$(basename $R1 | sed s/_R1.fq.gz//)

##trim
trim_galore --paired $R1 $R2

##map
bowtie2 -p 4 --no-unal -x /data/Blizard-BrancoLab/Genomes/GRCh38/hg38 -1 ${base_name}_R1_val_1.fq.gz -2 ${base_name}_R2_val_2.fq.gz -S ${base_name}.sam

##filter based on MAPQ
samtools view -b -q 2 $base_name.sam > $base_name.bam

##deduplicate
samtools fixmate -@ 3 -m $base_name.bam $base_name-fixmate.bam
samtools sort -@ 3 -o $base_name-sorted.bam $base_name-fixmate.bam
samtools markdup -@ 3 -r $base_name-sorted.bam $base_name-dedup.bam
samtools index -@ 3 $base_name-dedup.bam

##make bigwig
bamCoverage --bam $base_name-dedup.bam -o $base_name.bw --binSize 200 --normalizeUsing CPM
