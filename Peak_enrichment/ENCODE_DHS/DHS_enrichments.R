library(tidyverse)


##output from permutation pipeline

dhs = as_tibble(t(read.delim('DHS_enrichments.txt', row.names=1)), rownames='File.accession')


##ENCODE metadata

meta = read_tsv('metadata.tsv') %>% filter(grepl('released', File.Status))
dhs2 = inner_join(dhs, meta, by='File.accession')


##plot enrichments for a given family (Figure 1C and Supplementary Figure 1B)

plot.dhs = function(te) {
	quartz(w=2.5,h=2.5)
	select(dhs2,tissue=Biosample.term.name,enr=te) %>%
	ggplot(aes(x=tissue, y=enr, colour=tissue)) +
		#stat_summary(geom='bar',fun='mean', fill='wheat', colour='black', width=0.8) +
		geom_point(position=position_jitter(0.2), size=1, show.legend=FALSE) +
		theme_classic() +
		xlab('') +
		ylab('observed/expected')
}

plot.dhs('LTR10A')
plot.dhs('LTR3A')
plot.dhs('MER21A')
plot.dhs('MER41B')
plot.dhs('MER61E')
plot.dhs('LTR13')
plot.dhs('LTR4')


##select placenta-enriched families:
# 1. >80% of samples with enrichment above 1
# 2. minimum difference between median in placenta and other samples > 1

enriched = logical(ncol(dhs)-1)
names(enriched) = colnames(dhs)[-1]
for (i in 2:ncol(dhs)) {
	te = colnames(dhs)[i]
	enr = group_by(dhs2, Biosample.term.name) %>%
		summarise(p50 = sum(get(te)>1)/length(get(te)), med=median(get(te)))
	plac = enr$Biosample.term.name=='placenta'
	diff = enr$med[plac]-enr$med[!plac]
	if (min(diff>1) & enr$p50[plac]>0.8) enriched[i-1] = TRUE
}
