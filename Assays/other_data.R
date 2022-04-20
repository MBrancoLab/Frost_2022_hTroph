library(pzfx)
library(tidyverse)


##growth curve (Figure 6C)

df = read_pzfx('ENG.pzfx','ENG growth curves')
grow = tibble(sample = rep(gsub('_[[:digit:]]+$', '', colnames(df)[-1]), each=nrow(df)),
		day = rep(df$Var.1, ncol(df)-1),
		count = as.numeric(as.matrix(df[,-1]))) %>%
		filter(!is.na(count))

quartz(w=4.5,h=3)
ggplot(grow, aes(x=day, y=count, colour=sample)) +
	geom_line(stat='summary', fun='mean') +
	geom_point() +
	theme_classic() +
	xlab('Days in culture') +
	ylab('Cell count') +
	scale_colour_manual(values=c('#C00000','#F8CBAD'))

summary(aov(grow$count ~ grow$sample + grow$day))


##HLA-G FACS (Figure 6D)

df = read_pzfx('ENG.pzfx','HLAG+ cells following EvT differentiation')
df = df[,!is.na(df[1,])]
hlag = tibble(sample = gsub('_[[:digit:]]+$', '', colnames(df)),
	facs = as.numeric(df[1,]))

quartz(w=2.4, h=3)	
ggplot(hlag, aes(x=factor(sample, levels=unique(sample)), y=facs)) +
	geom_bar(stat='summary', fun=mean, col='black', fill='grey', width=0.6) +
	geom_point(position=position_jitter(0.1)) +
	xlab('') +
	ylab('% HLA-G+ cells') +
	theme_classic()
	

##sENG secretion (Figure 6E)

df = read_pzfx('ENG.pzfx','Soluble Endoglin protein secreted by Syncytiotrophoblast')
df = df[,!is.na(df[1,])]
sec = tibble(sample = gsub('_[[:digit:]]+$', '', colnames(df)),
	seng = as.numeric(df[1,]))

quartz(w=2.4, h=3)	
ggplot(sec, aes(x=factor(sample, levels=unique(sample)), y=seng)) +
	geom_bar(stat='summary', fun=mean, col='black', fill='grey', width=0.6) +
	geom_point(position=position_jitter(0.1)) +
	xlab('') +
	ylab('[sENG] pg/ml') +
	theme_classic()

wilcox.test(sec$seng ~ sec$sample)

