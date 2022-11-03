library(tidyverse)


##read FIMO files

data = tibble()
for (f in list.files('./FIMO_inactive', full.names=TRUE)) {
	data = rbind(data, read_tsv(f, comment='#'))
}


##filter for q<0.2

sel.data = data %>% filter(`q-value`<0.2) %>%
	mutate(te=gsub('_[[:digit:]]+$', '', sequence_name))


##summarise

data.sum = sel.data %>% group_by(te, motif_alt_id) %>%
	summarise(count=length(unique(sequence_name)))


##count number of inactive TEs per family

te = scan('../Active_families/te_families.txt', character())
total = numeric(length(te))

for (i in 1:length(te)) {
	system(paste('~/Documents/bedtools2/bin/intersectBed -v -a ../Annotations/',te[i],'_hg38.bed -b ../Peaks/hTSC_CutTag_H3K27ac.relaxed.bed > temp.bed', sep=''))
	total[i] = nrow(read.delim('temp.bed', header=F))
	unlink('temp.bed')
}
inact.sum = tibble(te, total)


##calculate percentage

data.sum = inner_join(data.sum, inact.sum, by='te') %>% mutate(perc=count/total*100)


##plot full table (Supplementary Figure 2B)

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

