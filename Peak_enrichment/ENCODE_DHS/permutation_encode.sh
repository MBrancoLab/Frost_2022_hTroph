#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_rt=12:0:0
#$ -l h_vmem=2G
#$ -t 1-133

module load bedtools/2.26.0
module load R

##download DHS file
url=$(sed -n "${SGE_TASK_ID}p" bed_list.txt)
wget $url

##run permutation pipeline
file=$(basename $url)
Rscript permutation_encode.R $file
