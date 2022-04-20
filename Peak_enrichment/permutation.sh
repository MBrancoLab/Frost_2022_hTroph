#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_rt=12:0:0
#$ -l h_vmem=2G
#$ -t 1-13

module load bedtools/2.26.0
module load R

##get peak file
file=$(sed -n "${SGE_TASK_ID}p" peak_list.txt)

##run permutation pipeline
Rscript permutation_hg38.R $file
