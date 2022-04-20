#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_rt=1:0:0
#$ -l h_vmem=2G
#$ -t 1-21

module load python/3.6.3
module load use.own
module load localtools

te=$(sed -n "${SGE_TASK_ID}p" te_families.txt)

computeMatrix reference-point -S hTSC_CutTag_cJun.bw hTSC_CutTag_JunD.bw \
	-R ${te}_hg38.bed \
	--referencePoint "center" \
	-a 3000 \
	-b 3000 \
	-bs 50 \
	-o $te.mat

plotProfile --matrixFile $te.mat \
	--outFileName $te.svg \
	--outFileNameData $te-profile.txt \
	--averageType "mean"
	--plotHeight 3 \
	--plotWidth 3
