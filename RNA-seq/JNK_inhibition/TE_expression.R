library(tidyverse)


##read subfamily counts

subf.files = list.files(pattern='subFcounts.txt')
subf = tibble()
for (f in subf.files) {
	data = read_tsv(f) %>% mutate(group = substring(Sample, 1, nchar(Sample)-1))
	subf = rbind(subf, data)
}

subf$te = unlist(lapply(strsplit(subf$`Subfamily:Family:Class`, split=':'), function(x) x[1]))


##normalise by 75th quantile

fac = group_by(subf, Sample) %>% summarise(q75 = quantile(tot_counts, 0.75))
subf = full_join(subf, fac, by='Sample') %>%
	mutate(norm = tot_counts/q75)


##plot selected subfamilies

fam = c(	'HERVIP10F-int','HERVIP10FH-int','Harlequin-int','HERVK3-int')

filter(subf, te %in% fam) %>%
ggplot(aes(x=factor(te, levels=fam), y=norm, fill=group)) +
	geom_bar(stat='summary', fun='mean', position=position_dodge()) +
	geom_point(position=position_dodge(0.9)) +
	theme_classic() +
	xlab('') +
	ylab('Normalised expression')

