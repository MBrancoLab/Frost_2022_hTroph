library(tidyverse)


##read data from permutation pipeline

pfiles = list.files(path='../Peak_enrichment', pattern='_overlap.txt')
sample = unlist(lapply(strsplit(pfiles,split='\\.'), function(x) x[1]))

over = tibble()
for (i in 1:length(pfiles)) {
	sub = read_tsv(paste('../Peak_enrichment', pfiles[i], sep='/')) %>% add_column(sample=sample[i])
	over = rbind(over, sub)
}
colnames(over)[1] = 'rep'


##highlight hTSC- and hESC-active TE families

ts.fam = scan('../Active_families/te_families.txt', character())
es.fam = scan('../Active_families/esc_families.txt', character())

over = mutate(over, activity = case_when(rep %in% ts.fam ~ 'TS',	
							rep %in% es.fam ~ 'ES',
							TRUE ~ '-'))


##enrichment scatter plots (Figure 2C and Supplementary Figure 2B)

enr.plot = function(exp) {	
	filter(over, sample==exp) %>% arrange(activity) %>%
	ggplot(aes(x=log10(real.overlap), y=real.overlap/random.mean)) +
		geom_point(aes(colour=activity, size=activity)) +
		theme_classic() +
		scale_color_manual(values=c('grey','blue','red')) +
		scale_size_manual(values=c(0.5,1,1)) +
		xlab('log10 n overlapped') +
		ylab('observed/expected')
}

quartz(w=4, h=3)
enr.plot('hTSC_CutTag_cJun') + ylim(0,40)
enr.plot('CT29_hTSC_GATA3') + ylim(0,20)
enr.plot('CT29_hTSC_TEAD4') + ylim(0,10)
enr.plot('CT29_hTSC_TFAP2C') + ylim(0,15)
enr.plot('hTSC_CutTag_JunD') + ylim(0,50)