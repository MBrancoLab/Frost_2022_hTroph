library(tidyverse)
library(circlize)


##intersectBed path

intersect.bed = '~/Documents/bedtools2/bin/intersectBed'


##intersect function

intersect.te = function(te) {
	
	#read TE annotation
	te.path = paste('../Annotations/', te, '_hg38.bed', sep='')
	bed = read_tsv(te.path,	col_names=c('chr','start','end','id','score','strand'))

	#intersect with H3K27ac	
	system(paste(intersect.bed, '-c -a', te.path, '-b ../Peaks/hTSC_CutTag_H3K27ac.relaxed.bed > temp.bed'))
	over = read.delim('temp.bed', header=F)
	bed = add_column(bed, H3K27ac = over$V7)

	#intersect with H3K4me1
	system(paste(intersect.bed, '-c -a', te.path, '-b ../Peaks/hTSC_CutTag_H3K4me1.relaxed.bed > temp.bed'))
	over = read.delim('temp.bed', header=F)
	bed = add_column(bed, H3K4me1 = over$V7)

	#intersect with H3K4me3
	system(paste(intersect.bed, '-c -a', te.path, '-b ../Peaks/hTSC_CutTag_H3K4me3.relaxed.bed > temp.bed'))
	over = read.delim('temp.bed', header=F)
	bed = add_column(bed, H3K4me3 = over$V7)

	#intersect with H3K9me3
	system(paste(intersect.bed, '-c -a', te.path, '-b ../Peaks/CT29_hTSC_H3K9me3.relaxed.bed > temp.bed'))
	over = read.delim('temp.bed', header=F)
	bed = add_column(bed, H3K9me3 = over$V7)

	#intersect with H3K27me3
	system(paste(intersect.bed, '-c -a', te.path, '-b ../Peaks/CT29_hTSC_H3K27me3.relaxed.bed > temp.bed'))
	over = read.delim('temp.bed', header=F)
	bed = add_column(bed, H3K27me3 = over$V7)
	unlink('temp.bed')
	
	return(bed)
}


##circos plotting function

circos.te = function(bed) {
	
	#order
	bed = arrange(bed, H3K27ac, H3K4me1, H3K4me3, H3K9me3, H3K27me3)

	#make plot
	circos.clear()
	circos.par(gap.degree=0, track.height=0.15)	
	circos.heatmap(bed$H3K27ac,col=colorRamp2(c(0, 1), c("white", "grey")), cluster=FALSE)
	circos.heatmap(bed$H3K4me1,col=colorRamp2(c(0, 1), c("white", "orange")))
	circos.heatmap(bed$H3K4me3,col=colorRamp2(c(0, 1), c("white", "green3")))
	circos.heatmap(bed$H3K9me3,col=colorRamp2(c(0, 1), c("white", "red1")))
	circos.heatmap(bed$H3K27me3,col=colorRamp2(c(0, 1), c("white", "purple")))
}


##get intersections

te.fams = c('LTR10A','LTR10F','LTR23','LTR2B','LTR3A','LTR7C','LTR8','LTR8B',
	'MER11D','MER21A','MER31A','MER39','MER39B','MER41A','MER41B','MER41C',
	'MER61D','MER61E')

intersect.list = list()
for (i in 1:length(te.fams)) {
	intersect.list[[i]] = intersect.te(te.fams[i])
}
names(intersect.list) = te.fams


##proportions of selected signatures

k27ac.k4me1 = unlist(lapply(intersect.list, function(x) {
		sum(x$H3K27ac>0 & x$H3K4me1>0 & x$H3K4me3==0)/sum(x$H3K27ac>0)
	}))

k27ac.k4me3 = unlist(lapply(intersect.list, function(x) {
		sum(x$H3K27ac>0 & x$H3K4me1==0 & x$H3K4me3>0)/sum(x$H3K27ac>0)
	}))

k4me1 = unlist(lapply(intersect.list, function(x) {
		sum(x$H3K27ac==0 & x$H3K4me1>0 & x$H3K4me3==0)/sum(x$H3K4me1>0)
	}))


##plots (Figure 1E)

quartz(w=3,h=3)
circos.te(intersect.list$LTR10A)
circos.te(intersect.list$LTR3A)
circos.te(intersect.list$MER21A)
circos.te(intersect.list$MER41B)
circos.te(intersect.list$MER61C)
