#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_rt=1:0:0
#$ -l h_vmem=2G

module load python/3.6.3
module load use.own
module load localtools


#Figure 1C and Supplementary Figure 1B

while read te; do

computeMatrix reference-point -S hTSC_CutTag_H3K27ac-2.bw SRR8383480_H1_ESC_H3K27ac-1.bw \
	-R ${te}_hg38.bed \
	--referencePoint "center" \
	-a 3000 \
	-b 3000 \
	-bs 50 \
	-o $te.mat

plotHeatmap --matrixFile $te.mat \
	--outFileName $te.svg \
	--sortUsing sum \
	--heatmapHeight 12 \
	--heatmapWidth 3 \
	--whatToShow "heatmap and colorbar"

done < te_families.txt
