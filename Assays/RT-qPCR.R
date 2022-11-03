library(pzfx)
library(tidyverse)


##functions to parse prism data

pzfx_to_tibble = function(file, table) {
	df = read_pzfx(file, table)
	df = df[,!is.na(df[1,])]
	tb = tibble(sample = gsub('_[[:digit:]]+$', '', colnames(df)),
		expr = as.numeric(df[1,]))
	return(tb)
}

split_sample = function(.data) {
	mutate(.data, cells = unlist(lapply(strsplit(sample, split=' '), function(x) x[1])),
		sgrna = unlist(lapply(strsplit(sample, split=' '), function(x) x[2])))
}


##basic plot functions

plot_single = function(data) {
	ggplot(data, aes(x=factor(sample, levels=unique(sample)), y=expr)) +
		geom_bar(stat='summary', fun=mean, col='black', fill='lightblue', width=0.6) +
		geom_point(position=position_jitter(0.1)) +
		xlab('') +
		ylab('Relative expression') +
		theme_classic()
}

plot_multi = function(data) {
	ggplot(data, aes(x=factor(cells, levels=unique(cells)), y=expr, fill=factor(sgrna, levels=unique(sgrna)))) +
		geom_bar(stat='summary', fun=mean, col='black', width=0.6, position=position_dodge(), show.legend=FALSE) +
		geom_point(position=position_jitterdodge(jitter.width=0.1,dodge.width=0.6), show.legend=FALSE) +
		xlab('') +
		ylab('Relative expression') +
		theme_classic()
}

plot_diff = function(file, table) {
	data = pzfx_to_tibble(file,table) %>%
	mutate(cells = unlist(lapply(strsplit(sample, split=' '), function(x) x[3])),
		sgrna = unlist(lapply(strsplit(sample, split=' '), function(x) x[1])))
	plot_multi(data)
}


##ADAM9 excision (Figure 5A and Supplementary Figure 5A)

pzfx_tables('ADAM9.pzfx')

adam9.p1 = pzfx_to_tibble('ADAM9.pzfx', 'XY: ADAM9_iso12')
adam9.p2 = pzfx_to_tibble('ADAM9.pzfx', 'XY: ADAM9_iso23')
htra4 = pzfx_to_tibble('ADAM9.pzfx', 'XY: HTRA4')
tm2d2 = pzfx_to_tibble('ADAM9.pzfx', 'XY: TM2D2')

quartz(w=2.4, h=3)
plot_single(adam9.p1)
wilcox.test(adam9.p1$expr ~ adam9.p1$sample)

quartz(w=2.4, h=3)
plot_single(adam9.p2)
wilcox.test(adam9.p2$expr ~ adam9.p2$sample)

quartz(w=2.4, h=3)
plot_single(htra4)
wilcox.test(htra4$expr ~ htra4$sample)

quartz(w=2.4, h=3)
plot_single(tm2d2)
wilcox.test(tm2d2$expr ~ tm2d2$sample)


##CSF1R excision (Figure 5B and Supplementary Figure 5C)

pzfx_tables('CSF1R.pzfx')

csf1r = rbind(pzfx_to_tibble('CSF1R.pzfx', 'LTR10A at CSF1R hTSC all+placenta BTS5'),
	pzfx_to_tibble('CSF1R.pzfx','LTR10A at CSF1R EVT all+placenta BTS5'),
	pzfx_to_tibble('CSF1R.pzfx','LTR10A at CSF1R both all+placenta - ST(3D) BTS5')) %>%
	split_sample()

csf.all = filter(csf1r, grepl('all', sample))
quartz(w=4, h=3)
plot_multi(csf.all)
wt.all = group_by(csf.all, cells) %>% summarise(p=wilcox.test(expr~sgrna)$p.value)
p.adjust(wt.all$p, method='BH')

csf.plac = filter(csf1r, grepl('placenta', sample))
quartz(w=4, h=3)
plot_multi(csf.plac)
wt.plac = group_by(csf.plac, cells) %>% summarise(p=wilcox.test(expr~sgrna)$p.value)
p.adjust(wt.plac$p, method='BH')

quartz(w=4, h=3)
plot_diff('CSF1R.pzfx', 'CSF1R TEAD4')
quartz(w=2.7, h=3)
plot_diff('CSF1R.pzfx', 'CSF1R SDC')
plot_diff('CSF1R.pzfx', 'CSF1R hCGb')
plot_diff('CSF1R.pzfx', 'CSF1R MMP')
plot_diff('CSF1R.pzfx', 'CSF1R HLAG')


##TWIST excision (Figure 5C and Supplementary Figure 5B)

pzfx_tables('TWIST.pzfx')

twist = pzfx_to_tibble('TWIST.pzfx','TWIST isoform 2 differentiated')
twist$sample[twist$sample=='MER41B excision'] = 'hTSC MER41B excision'
twist = twist %>% split_sample()

twist.tsc = filter(twist, cells=='hTSC')
quartz(w=2.4, h=3)
plot_single(twist.tsc)

quartz(w=4, h=3)
plot_multi(twist)
wt.twist = group_by(twist, cells) %>% summarise(p=wilcox.test(expr~sgrna)$p.value)
p.adjust(wt.twist$p, method='BH')


##PSG5 excision (Figure 5D and Supplementary Figure 5D)

pzfx_tables('PSG5.pzfx')

psg5 = pzfx_to_tibble('PSG5.pzfx','XY: PSG5_iso36') %>% split_sample()
psg5 = psg5[order(match(psg5$cells, c('hTSC','EvT','ST(3D)'))),]

quartz(w=4, h=3)
plot_multi(psg5)
wt.psg = group_by(psg5, cells) %>% summarise(p=wilcox.test(expr~sgrna)$p.value)
p.adjust(wt.psg$p, method='BH')

quartz(w=4, h=3)
plot_diff('PSG5.pzfx', 'TEAD4_Branco')
quartz(w=2.7, h=3)
plot_diff('PSG5.pzfx', 'SDC expression ST(3D) differentiation')
plot_diff('PSG5.pzfx', 'hCGb expression ST(3D) differentiation')
plot_diff('PSG5.pzfx', 'MMP2 EvT_Branco')
plot_diff('PSG5.pzfx', 'HLAG EvT_Branco')


##ENG excision (Figure 6B,D,E and Supplementary Figure 6A,C)

pzfx_tables('ENG.pzfx')

eng.p1 = pzfx_to_tibble('ENG.pzfx','Eng new non normalised common')
quartz(w=3, h=4)
plot_single(eng.p1)
p = c(with(filter(eng.p1,sample!='Non-LTR excision'), wilcox.test(expr~sample)$p.value),
	with(filter(eng.p1,sample!='LTR10A excision'), wilcox.test(expr~sample))$p.value)
p.adjust(p, method='BH')

eng.p2 = pzfx_to_tibble('ENG.pzfx','Eng new non normalised var1&2')
quartz(w=3, h=4)
plot_single(eng.p2)
p = c(with(filter(eng.p2,sample!='Non-LTR excision'), wilcox.test(expr~sample)$p.value),
	with(filter(eng.p2,sample!='LTR10A excision'), wilcox.test(expr~sample))$p.value)
p.adjust(p, method='BH')

eng.st = pzfx_to_tibble('ENG.pzfx','Eng new non normalised var1&2 ST(3D)')
quartz(w=2.4, h=3)
plot_single(eng.st)
wilcox.test(eng.st$expr ~ eng.st$sample)

eng.evt = pzfx_to_tibble('ENG.pzfx','Eng new non normalised var1&2 EvT')
quartz(w=2.4, h=3)
plot_single(eng.evt)
wilcox.test(eng.evt$expr ~ eng.evt$sample)

ak1 = pzfx_to_tibble('ENG.pzfx','AK1_ENG_FLANK')
fpgs = pzfx_to_tibble('ENG.pzfx','FPGS_ENG_FLANK')
quartz(w=3, h=3)
plot_single(ak1)
plot_single(fpgs)

quartz(w=2.7, h=3)
plot_diff('ENG.pzfx', 'ENG_diff_TEAD4')
plot_diff('ENG.pzfx', 'ENG_diff_SDC')
plot_diff('ENG.pzfx', 'ENG_diff_CSH')


##SR11302 treatment (Supplementary Figure 4??)

sr = read_tsv('SR11302.txt')

quartz(w=2.4, h=3)
plot_single(filter(sr, gene=='MMP14'))
plot_single(filter(sr, gene=='NOS3'))
plot_single(filter(sr, gene=='ENG'))


