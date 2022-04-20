#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=2G
#$ -t 1-9

module load use.own
module load seacr
module load bedtools
module load R

ip=$(sed -n "${SGE_TASK_ID}p" seacr_ip.txt)
input=$(sed -n "${SGE_TASK_ID}p" seacr_input.txt)
base=$(echo $ip | sed s/.fragments.bedgraph//)

bash SEACR_1.2.sh $ip $input norm relaxed $base
