library(tidyverse)


##read TE transcript data

flist = c('Control1_L7.D708', 'Control2_L7.D709', 'Control3_L5.D710', 'Control4_L4.D711')
data = tibble()
for (f in flist) {
	data = rbind(data,
		read_tsv(paste('te_promoters/', f, '-TEtranscripts.txt', sep=''), col_names=FALSE) %>%
			mutate(te = paste(X1,X2,X3,X4, sep='_')) %>%
			add_column(sample = f)
		)
}


##filter for selected families

fam = scan('../Active_families/te_families.txt', character())
sel = filter(data, X4 %in% fam)


##remove redundancy

nonred = sel %>% group_by(te, X10, X12, sample) %>% summarise(fpkm = max(X11)) %>%
	pivot_wider(names_from = sample, values_from = fpkm)

colnames(nonred)[2:3] = c('gene', 'strand')


##write

write_tsv(nonred, 'TE_promoters.txt')