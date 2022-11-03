library(tidyverse)


##read data from permutation pipeline

pfiles = list.files(pattern='_overlap.txt')
sample = unlist(lapply(strsplit(pfiles,split='\\.'), function(x) x[1]))

over = tibble()
for (i in 1:length(pfiles)) {
	sub = read_tsv(pfiles[i]) %>% add_column(sample=sample[i])
	over = rbind(over, sub)
}
colnames(over)[1] = 'rep'


##get enriched repeats

enriched = mutate(over, is.enr = pval<0.05 & log.enrich>=1 & real.overlap>=10) %>%
	select(rep, sample, is.enr) %>%
	pivot_wider(names_from=sample, values_from=is.enr)


##select H3K27ac-enriched repeats

k27.reps = filter(enriched, hTSC_CutTag_H3K27ac, hTSC_ChIP_H3K27ac_peaks, Cytotrophoblast_H3K27ac_peaks) %>%
	pull(rep)
over.k27 = mutate(over, is.k27 = rep %in% k27.reps)

#for ES-enriched repeats:
#filter(enriched, H1_ESC_H3K27ac, !(hTSC_CutTag_H3K27ac | hTSC_ChIP_H3K27ac_peaks | Cytotrophoblast_H3K27ac_peaks)) %>% pull(rep)


##enrichment scatter plots (Figure 1B and Supplementary Figure 1A)

enr.scatter = function(sample_to_plot, ymax=15) {
	filter(over.k27, sample==sample_to_plot) %>% arrange(is.k27) %>%
	ggplot(aes(x=log10(real.overlap), y=real.overlap/random.mean)) +
		geom_point(aes(colour=is.k27, size=is.k27)) +
		theme_classic() +
		scale_color_manual(values=c('grey','black')) +
		scale_size_manual(values=c(0.5,1)) +
		ylim(0,ymax) +
		xlab('log10 n overlapped') +
		ylab('observed/expected')
}

quartz(w=4, h=3)
enr.scatter('hTSC_CutTag_H3K27ac')
enr.scatter('H1_ESC_H3K27ac')
enr.scatter('hTSC_ChIP_H3K27ac_peaks', ymax=32)
enr.scatter('Cytotrophoblast_H3K27ac_peaks')


##write K27ac-enriched TE family lists into file

out.table = filter(over, sample=='hTSC_CutTag_H3K27ac', rep %in% k27.reps) %>% select(rep,total)

out.samples = c('hTSC_CutTag_H3K27ac', 'hTSC_CutTag_H3K4me1', 'hTSC_CutTag_H3K4me3',
	'hTSC_ChIP_H3K27ac_peaks', 'hTSC_ChIP_H3K4me3_peaks',
	'Cytotrophoblast_H3K27ac_peaks', 'Cytotrophoblast_H3K4me3_peaks',
	'CT29_hTSC_H3K9me3', 'CT27_hTSC_H3K27me3')
for (s in out.samples) {
	out.table = full_join(out.table,
		filter(over, sample==s, rep %in% k27.reps) %>%
		mutate(enrichment = real.overlap/random.mean) %>%
		select(rep, enrichment, pval,real.overlap),
		by='rep')
}

sample.header = paste(out.samples,collapse='\t\t\t')
write(paste('\t\t',sample.header), 'k27ac_families.txt', sep='')
write.table(out.table, 'k27ac_families.txt', append=TRUE,
	sep='\t', quote=FALSE, row.names=FALSE)

