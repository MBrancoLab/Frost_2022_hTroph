library(tidyverse)


##read data from permutation pipeline

over = rbind(read_tsv('../Peak_enrichment/hTSC_CutTag_H3K27ac.relaxed.bed_overlap.txt') %>%
					add_column(cell='hTSC'),
			read_tsv('../Peak_enrichment/BTS5_EVT_H3K27ac.relaxed.bed_overlap.txt') %>%
					add_column(cell='EVT'))
colnames(over)[1] = 'rep'


##enrichment change in selected families (Supplementary Figure 1E)

te.fams = c('LTR10A','LTR10F','LTR23','LTR2B','LTR3A','LTR7C','LTR8','LTR8B',
	'MER11D','MER21A','MER31A','MER39','MER39B','MER41A','MER41B','MER41C',
	'MER61D','MER61E')
	
sel = filter(over, rep %in% te.fams) %>%
	select(rep, log.enrich, cell) %>%
	pivot_wider(values_from=log.enrich, names_from=cell)

ggplot(sel, aes(x=hTSC, y=rep)) +
	geom_vline(xintercept=1, linetype='dashed', col='grey') +
	geom_segment(aes(x=hTSC, y=rep, xend=EVT, yend=rep))+
	geom_point(size=2, col='blue') +
	geom_point(aes(x=EVT, y=rep), size=2, col='red') +
	theme_classic() +
	xlab('log2 obs/exp') +
	ylab('')


##concatenate EVT-enriched TEs

evt.tes = c('LTR23','LTR2B','LTR3A','LTR8','MER11D','MER21A','MER31A','MER41B','MER41C')
evt.files = paste(paste('../Annotations/', evt.tes, '_hg38.bed', sep=''),collapse=' ')
system(paste('cat', evt.files, '> evt_tes.bed'))
bed = read_tsv('evt_tes.bed', col_names=c('chr','start','end','id','score','strand'))


##intersect with hTSC and EVT peaks

intersect.bed = '~/Documents/bedtools2/bin/intersectBed'

system(paste(intersect.bed, '-c -a evt_tes.bed -b ../Peaks/hTSC_CutTag_H3K27ac.relaxed.bed > temp.bed'))
temp = read.delim('temp.bed', header=F)
bed = add_column(bed, H3K27ac = temp$V7>0)

system(paste(intersect.bed, '-c -a evt_tes.bed -b ../Peaks/hTSC_CutTag_H3K4me1.relaxed.bed > temp.bed'))
temp = read.delim('temp.bed', header=F)
bed = add_column(bed, H3K4me1 = temp$V7>0)

system(paste(intersect.bed, '-c -a evt_tes.bed -b ../Peaks/hTSC_CutTag_H3K4me3.relaxed.bed > temp.bed'))
temp = read.delim('temp.bed', header=F)
bed = add_column(bed, H3K4me3 = temp$V7>0)

system(paste(intersect.bed, '-c -a evt_tes.bed -b ../Peaks/BTS5_EVT_H3K27ac.relaxed.bed > temp.bed'))
temp = read.delim('temp.bed', header=F)
bed = add_column(bed, EVT = temp$V7>0)

unlink(c('evt_tes.bed', 'temp.bed'))


##classify hTSC states

bed = mutate(bed, TSC.state = case_when(H3K27ac & H3K4me1 & !H3K4me3 ~ 'active_enh',
			!H3K27ac & H3K4me1 & !H3K4me3 ~ 'poised_enh',
			H3K4me3 ~ 'promoter',
			!H3K27ac & !H3K4me1 & !H3K4me3 ~ 'inactive'),
		repname = gsub('_[[:digit:]]+','',id))


##get frequencies

freq = group_by(bed, repname, TSC.state, EVT) %>% summarise(frequency = length(EVT))


##plot TSC state frequencies for EVT+ TEs (Supplementary Figure 1F)

filter(freq, EVT, TSC.state!='NA') %>%
ggplot(aes(x=frequency, y=repname, fill=TSC.state)) +
	geom_bar(stat='identity') +
	theme_classic() +
	ylab('') +
	xlab('Number of EVT-active TEs')





