#!/bin/sh
#$ -cwd
#$ -pe smp 4
#$ -l h_rt=12:0:0
#$ -l h_vmem=12G
#$ -t 1-6


r1=$(sed -n "${SGE_TASK_ID}p" fastq_list.txt)
r2=$(echo $r1 | sed s/R1.fastq.gz/R2.fastq.gz/)
base=$(echo $r1 | sed s/-R1.fastq.gz//)


module load anaconda2

rlength=$(zcat $r1 | head | sed -n "2p" | tr -d '\n' | wc -m)

##activate virtual environment
source activate squire

##map
squire Map -1 $r1 -2 $r2 -f ~/genomes/GRCh38/squire_fetch -r $rlength -n $base -o ${base}_map -p 4 -v

##count
squire Count -m squire_map -c ~/genomes/GRCh38/squire_clean -f ~/genomes/GRCh38/squire_fetch -r $rlength -n $base -m ${base}_map -o ${base}_count -p 4 -s 1 -v


##draw
squire Draw -f ~/genomes/GRCh38/squire_fetch -m ${base}_map -o ${base}_draw -s 1 -b hg38
source deactivate squire

