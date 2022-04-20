library(tidyverse)
source('get_nearest_TE.R')


##Okae RNA data

rna = read.delim('RNA_hTSC_hg38.txt.gz', as.is=TRUE)


##Re-normalise by 75th quantile
#(particularly relevant for ST, which have massively outlier genes)

fpkm = 2^rna[,7:ncol(rna)] - 1

q75 = apply(fpkm, 2, function(x) quantile(x, 0.75))
cf = q75/mean(q75)
norm = fpkm/matrix(rep(cf, nrow(fpkm)), nrow=nrow(fpkm), byrow=T)

rna[,7:ncol(rna)] = log2(norm+1)


##make tibble with average expression

rna.av = tibble(gene = rep(rna$gene,7),
				id = rep(rna$id, 7),
				chr = rep(rna$chr, 7),
				start = rep(rna$start, 7),
				end = rep(rna$end, 7),
				strand = rep(rna$strand, 7),
				expr=c(rowMeans(rna[,grep('Stroma',colnames(rna))]),
					   rowMeans(rna[,grep('^CT',colnames(rna))]),
					   rowMeans(rna[,grep('^TS',colnames(rna))]),
					   rowMeans(rna[,grep('^EVT\\.\\.',colnames(rna))]),
					   rowMeans(rna[,grep('^EVT.TS',colnames(rna))]),
					   rowMeans(rna[,grep('^ST\\.\\.',colnames(rna))]),
					   rowMeans(rna[,grep('^ST.3D',colnames(rna))])),
				cells=rep(c('Stroma','CT','TSC','EVT','TSC.EVT','ST','TSC.ST'), each=nrow(rna)))


##get nearest TEs

rna.te = get_nearest_TE(rna.av, '../Annotations/H3K27ac_TEs.bed')
rna.te$distance[rna.te$distance==-1] = max(rna.te$distance)

rna.es = get_nearest_TE(rna.av, '../Annotations/ESC_H3K27ac_TEs.bed')
rna.es$distance[rna.es$distance==-1] = max(rna.es$distance)

rna.te = mutate(rna.te, dist.es = rna.es$distance,
			dgroup.es = cut(dist.es, c(0,10000,50000,100000, max(dist.es)), include.lowest=T))


##get expression relative to stroma

rel.expression = function(data, ref.cells) {
	ref = filter(data, cells==ref.cells)
	rel = add_column(data, ref = ref$expr[match(data$id, ref$id)]) %>%
		mutate(rel=expr-ref, min.expr = expr>1 | ref>1) %>%
		filter(cells!=ref.cells)
	return(rel)
}

rna.rel = rel.expression(rna.te, 'Stroma') %>%
	mutate(dgroup = cut(distance, c(0,10000,50000,100000, max(distance)), include.lowest=T))


##plot relative expression as a function of distance (Figure 3A and Supplementary Figure 3B)

rel.plot = function(cell.selection, group.var='dgroup') {
	quartz(w=4.5,h=3)
	rna.rel %>% filter(cells %in% cell.selection, min.expr) %>%
	ggplot(aes(x=cells, y=rel, fill=get(group.var))) +
		geom_boxplot(outlier.colour=NA) +
		theme_classic() +
		ylim(-5,5) +
		scale_fill_manual(values=colorRampPalette(c('lightblue','darkblue'))(4)) +
		geom_hline(yintercept=0, linetype='dashed')
}

rel.plot(c('TSC','TSC.EVT','TSC.ST'))
rel.plot(c('CT','EVT','ST'))

rel.plot(c('TSC','TSC.EVT','TSC.ST'), group.var='dgroup.es')
rel.plot(c('CT','EVT','ST'), group.var='dgroup.es')


##ANOVA

rel.aov = function(cell.selection, group.var='dgroup') {
	aov.data = rna.rel %>% filter(cells %in% cell.selection, min.expr)
	model = aov(rel ~ cells * get(group.var), data=aov.data)
	all.p = TukeyHSD(model)[[3]]
	split.comp = strsplit(rownames(all.p), split='[-:]')
	same.cell = unlist(lapply(split.comp, function(x) x[1]==x[3]))
	return(all.p[same.cell,])
}

rel.aov(c('TSC','TSC.EVT','TSC.ST'))
rel.aov(c('CT','EVT','ST'))

rel.aov(c('TSC','TSC.EVT','TSC.ST'), group.var='dgroup.es')
rel.aov(c('CT','EVT','ST'), group.var='dgroup.es')



#################################
# ES-to-TS transdifferentiation #
#################################


##read Dong et al data

dong = read.delim('Dong_RNA_transDiff.txt.gz', as.is=T)


##normalise and make tibble

get.rpm = function(x) log2(x/sum(x)*1e6 + 0.1)

tran = tibble(gene = rep(dong$name, 2),
			cell = rep(c('H9','AN'), each=nrow(dong)),
			lfc = c(get.rpm(dong$H9_naive_TSC) - get.rpm(dong$H9_naive),
					get.rpm(dong$AN_naive_TSC) - get.rpm(dong$AN_naive)),
			ts.close = gene %in% rna.te$gene[rna.te$distance<50000],
			es.close = gene %in% rna.te$gene[rna.te$dist.es<50000])


##plot (Figure 3B)

to.plot = pivot_longer(tran, cols=c('ts.close','es.close'), names_to='close', values_to='is.close') %>%
	filter(is.close)
	
quartz(w=3, h=3)
ggplot(to.plot, aes(x=cell, y=lfc, fill=close)) +
	geom_boxplot(outlier.colour=NA) +
	theme_classic() +
	ylim(-5,5) +
	xlab('') +
	ylab('log2 FC')

TukeyHSD(aov(lfc ~ cell * close, data=to.plot))


############################
# Variably expressed genes #
############################

library(gplots)


##select genes with nearby TEs

te.expr = rna.te %>% filter(distance<50000) %>%
	pivot_wider(values_from=expr, names_from=cells) %>%
	select(gene, Stroma, TSC, TSC.EVT, TSC.ST)


##select variably expressed genes

te.var = te.expr %>% mutate(var = apply(select(te.expr, !gene), 1, var)) %>%
	filter(var>1)


##k-means clustering

te.mat = as.matrix(select(te.var, !c(gene, var)))
rownames(te.mat) = te.var$gene

set.seed(208)
k = kmeans(te.mat, centers=10)
k.order = c(1,5,10,3,6,8,2,4,7,9)


##heatmap (Figure 3C)

quartz(w=3.5, h=4)
heatmap.2(te.mat[order(match(k$cluster, k.order)),],
	col=hcl.colors, breaks=seq(0,9,0.2), 
	scale='none', Rowv=FALSE, Colv=FALSE, dendrogram='none',
	trace='none', density.info='none', cexCol=0.8)



#####################
# Specific families #
#####################

plot.fam = function(fam) {
	near.fam = get_nearest_TE(rna.av, '../Annotations/H3K27ac_TEs.bed', family=fam)
	near.fam$distance[near.fam$distance==-1] = max(near.fam$distance)
	rel.fam = rel.expression(near.fam, 'Stroma')
	sub = filter(rel.fam, distance<50000, min.expr, cells %in% c('TSC','TSC.EVT','TSC.ST'))
	
	quartz(w=2, h=3)
	gp = ggplot(sub, aes(x=cells, y=rel, color=cells)) +
		geom_boxplot(outlier.color=NA, show.legend=FALSE) +
		geom_point(size=0.8, position=position_jitter(0.2), show.legend=FALSE) +
		geom_hline(yintercept=0, linetype='dashed') +
		theme_classic() +
		xlab('') +
		ylab('Expression relative to stroma')
	print(gp)
	
	p = group_by(sub, cells) %>% summarise(p = wilcox.test(rel)$p.value)
	return(p.adjust(p$p, method='BH'))
}

#(Supplementary Figure 3C)
plot.fam('LTR8B')
plot.fam('MER11D')
plot.fam('LTR23')
plot.fam('MER41B')

