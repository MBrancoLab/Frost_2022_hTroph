#!/bin/sh
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=5:0:0
#$ -t 1-4

module load python/3.6.3
module load use.own
module load macs2


ip=$(sed -n "${SGE_TASK_ID}p" macs2_ip.txt)
input=$(sed -n "${SGE_TASK_ID}p" macs2_input.txt)
exp=$(echo $ip | awk -F- '{print $1}')

##call peaks (PE data, broad peaks)
macs2 callpeak -t $ip -c $input --broad -g hs -n $exp -q 0.05 -f BAMPE
