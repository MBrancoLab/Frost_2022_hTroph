library(tidyverse)


##read FIMO files

data = tibble()
for (f in list.files('./FIMO', full.names=TRUE)) {
	data = rbind(data, read_tsv(f, comment='#'))
}


##filter for selected motifs and q<0.2

sel.motif = read_tsv('selected_motifs.txt', col_names=c('id','name'))

sel.data = data %>% filter(motif_id %in% sel.motif$id, `q-value`<0.2) %>%
	mutate(te=unlist(lapply(strsplit(sequence_name,'-'), function(x) x[1])))


##summarise

data.sum = sel.data %>% group_by(te, motif_alt_id) %>%
	summarise(count=length(unique(sequence_name)))


##count number of active TEs per family

active = read_tsv('../Annotations/H3K27ac_TEs.bed',
	col_names=c('chr','start','end','te','score','strand'))

act.sum = active %>% group_by(te) %>% summarise(total=length(te))


##calculate percentage

data.sum = inner_join(data.sum, act.sum, by='te') %>% mutate(perc=count/total*100)


##plot full table (Supplementary Figure 2A)

quartz(w=10,h=4)
data.sum %>%
ggplot(aes(te,motif_alt_id,colour=te)) +
	theme_classic() + 
	xlab('') +
	ylab('') +
	geom_point(aes(size=perc), stroke=0) +
	scale_size_continuous(range = c(0.5,7)) +
	ylim(sort(unique(data.sum$motif_alt_id), decreasing=TRUE)) +
	guides(colour='none')


##plot selected TEs (Figure 2A)

motif.order = sort(unique(data.sum$motif_alt_id), decreasing=TRUE)

quartz(w=5.7, h=4)
data.sum %>% filter(te %in% c('LTR10A','LTR2B','LTR3A','LTR7C','LTR23','MER21A','MER41C','MER61E','MER11D')) %>%
ggplot(aes(te,motif_alt_id,colour=te)) +
	theme_classic() + 
	xlab('') +
	ylab('') +
	geom_point(aes(size=perc), stroke=0) +
	scale_size_continuous(range = c(0.5,7)) +
	ylim(motif.order[motif.order!='SRF']) +
	guides(colour='none')

