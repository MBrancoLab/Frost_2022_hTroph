library(tidyverse)
source('../RNA-seq/get_nearest_TE.R')


####################
# Get RNA-seq data #
####################


##function to parse Seqmonk output

parse_seqmonk = function(file) {
	df = read.delim(file)
	out = tibble(chr = paste('chr', df$Chromosome, sep=''),
				start = df$Start,
				end = df$End,
				strand = df$Feature.Strand,
				gene = df$Feature,
				id = df$ID,
				df[,13:ncol(df)])
	colnames(out)[7:ncol(out)] = colnames(df)[13:ncol(df)]
	return(out)
}


##Okae hTSC data

human = parse_seqmonk('hTSC_RNA_RPKM.txt.gz')
hs.te = get_nearest_TE(human, '../Annotations/H3K27ac_TEs.bed')
hs.te$distance[hs.te$distance==-1] = max(hs.te$distance)


##Cambuli mTSC data

mouse = parse_seqmonk('mTSC_RNA_RPKM.txt.gz')
mm.te = get_nearest_TE(mouse, 'mTSC_H3K27ac_TEs.bed')
mm.te$distance[mm.te$distance==-1] = max(mm.te$distance)


##Schmidt macTSC data

mac = parse_seqmonk('macTSC_RNA_RPKM.txt.gz')



#################
# Join datasets #
#################


##function to join datasets for 1-to-1 gene orthologues

ortho.genes = read.delim('orthologous_genes.txt.gz')

join.sp = function(hs.rna, sp.rna, hs.id, sp.id) {
	hs.unique = table(hs.id)==1
	sp.unique = table(sp.id)==1
	ortho1 = hs.id %in% names(hs.unique)[hs.unique] & sp.id %in% names(sp.unique)[sp.unique]
	hs.ortho = filter(hs.rna, id %in% hs.id[ortho1]) %>%
		mutate(id2 = sp.id[ortho1][match(id, hs.id[ortho1])])
	
	hs.sp = inner_join(hs.ortho, sp.rna, by=c('id2'='id'))
	colnames(hs.sp)[1:5] = c('chr', 'start', 'end', 'strand', 'gene')
	return(hs.sp)
}


##function to re-normalise data based on quantiles

norm.q = function(data, hs.val, sp.val) {
	norm = filter(data, !is.na(hs.val)) %>%
		add_column(sp.norm = hs.val[match(rank(sp.val, ties.method='random'), rank(hs.val, ties.method='random'))],
			fold.diff = hs.val-sp.norm)
	return(norm)
}


##join human and mouse data

hs.mm = join.sp(hs.te, mm.te, ortho.genes$Gene.stable.ID, ortho.genes$Mouse.gene.stable.ID) %>%
	mutate(close = case_when(distance.x<50000 & distance.y>50000 ~ 'human',
							distance.x>50000 & distance.y<50000 ~ 'mouse',
							distance.x<50000 & distance.y<50000 ~ 'both',
							TRUE ~ 'neither'))
hs.mm.q = norm.q(hs.mm, hs.mm$hTSC.bam, hs.mm$mTSC.bam)


##join human and macaque data

hs.mac = join.sp(hs.te, mac, ortho.genes$Gene.stable.ID, ortho.genes$Macaque.gene.stable.ID)
hs.mac.q = norm.q(hs.mac, hs.mac$hTSC.bam, hs.mac$macTSC.bam)


##add TE orthology info to human-macaque set

ortho.tes = read.delim('../Active_families/orthologues/rheMac10_orthologues.bed', header=F)
hs.mac.ortho = mutate(hs.mac.q, is.ortho = paste(te.chr,te.start) %in% paste(ortho.tes$V1,ortho.tes$V2))



#########
# Plots #
#########


##human-macaque boxplot (Figure 4C)

mac.sel = filter(hs.mac.ortho, distance<100000, hTSC.bam>-1 | sp.norm>-1)

quartz(w=2, h=3)
ggplot(mac.sel, aes(x=is.ortho, y=fold.diff)) +
	geom_boxplot(outlier.color=NA, fill='tomato') +
	theme_classic() +
	ylim(-5,5) +
	geom_hline(yintercept=0, linetype='dashed') +
	xlab('') +
	ylab('log2 human/macaque')

wilcox.test(mac.sel$fold.diff ~ mac.sel$is.ortho)


##human-mouse boxplot (Figure 4D)

mm.sel = filter(hs.mm.q, hTSC.bam>-1 | sp.norm>-1, close!='both')

quartz(w=2.5, h=3)
ggplot(mm.sel, aes(x=close, y=fold.diff)) +
	geom_boxplot(outlier.color=NA, fill='tomato') +
	theme_classic() +
	ylim(-5,5) +
	geom_hline(yintercept=0, linetype='dashed') +
	xlab('') +
	ylab('log2 human/mouse')

TukeyHSD(aov(mm.sel$fold.diff ~ mm.sel$close))


##human-macaque individual families

#sel.te = c('LTR2B', 'LTR23', 'MER39B')
#mac.fam = filter(hs.mac.ortho, distance<100000, hTSC.bam>-1 | sp.norm>-1,
#				te.name %in% sel.te,
#				!(te.name=='LTR2B' & is.ortho)) #remove group with single data point

#quartz(w=5, h=3)
#ggplot(mac.fam, aes(x=factor(te.name, levels=sel.te), y=fold.diff, colour=is.ortho)) +
#	geom_violin() +
#	geom_point(position=position_jitterdodge(jitter.width=0.2)) +
#	theme_classic() +
#	geom_hline(yintercept=0, linetype='dashed') +
#	xlab('') +
#	ylab('log2 human/macaque')

#p = group_by(mac.fam, te.name, is.ortho) %>% summarise(p = wilcox.test(fold.diff)$p.value)
#p.adjust(p$p, method='BH')

