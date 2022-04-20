library('DESeq2')
library('topGO')
library('org.Hs.eg.db')
library('tidyverse')


##read raw counts

reads = read.table('AP1_read_counts.txt', header = T, row.names = 1)
colData = DataFrame(sample=factor(rep(c('DMSO','SP600125'),each=3)))


##get DE genes

dds = DESeqDataSetFromMatrix(reads,colData,design=~sample)
dds = DESeq(dds)
res = results(dds,contrast=c('sample','SP600125','DMSO'))
up = subset(res, padj<0.05 & log2FoldChange>1)
down = subset(res, padj<0.05 & log2FoldChange<(-1))


##gene ontology analysis (Supplementary Figure 4A)

go.ann = function(de.genes) {
	genes = factor(as.integer(rownames(res) %in% rownames(de.genes)))
	names(genes) = rownames(res)
	go = new('topGOdata', ontology = 'BP',
		allGenes = genes, geneSelectionFun = function(x) x==0,
		annotationFun = annFUN.org, mapping = 'org.Hs.eg.db', ID = 'SYMBOL')
	return(go)
}

go.plot = function(go, nodes=8, xmax=20) {
	rt = runTest(go, algorithm='classic', statistic='fisher')
	gt = GenTable(go, classicFisher=rt, topNodes=nodes)

	ggplot(gt, aes(x=-log10(as.numeric(classicFisher)), y=Term)) +
		geom_bar(stat='identity', fill='grey', colour='black', width=0.7) +
		xlim(0, xmax) +
		theme_classic() +
		ylab('') +
		xlab('-log10 pval')
}

go.up = go.ann(up)
go.down = go.ann(down)

go.plot(go.up)
go.plot(go.down)


##generate VST expression table

vsd = varianceStabilizingTransformation(dds)
expr = assay(vsd)
#write.table(expr,'AP1_RNA_vsd.txt',sep='\t',quote=F,col.names=NA) #for GEO submission


##Check expression of specific genes (Supplementary Figure 4B)

expr.tib = as_tibble(expr, rownames='gene') %>%
	pivot_longer(cols=!gene, values_to='log.expr', names_to='sample') %>%
	mutate(group=substr(sample,1,nchar(sample)-2), expr=2^log.expr)

plot.gene = function(gene.name) {
	filter(expr.tib, gene==gene.name) %>%
	ggplot(aes(x=group, y=expr)) +
		geom_bar(stat='summary', fun=mean, width=0.7, fill='lightblue', col='black') +
		geom_point(position=position_jitter(0.1)) +
		theme_classic() +
		xlab('') +
		ylab('Normalised expression')
}

quartz(w=2, h=3)
plot.gene('JUN')
plot.gene('MMP14')
plot.gene('NOS3')


##average expression, add gene coordinates

annot = read_tsv('../RNA_hTSC_hg38.txt.gz') %>% select(gene, id, chr, start, end, strand)
dmso = rowMeans(expr[,grep('DMSO',colnames(expr))])
sp = rowMeans(expr[,grep('SP600125',colnames(expr))])
gene.match = match(annot$gene, rownames(expr))
rna.av = add_column(annot, dmso = dmso[gene.match], sp = sp[gene.match]) %>%
	mutate(fc = sp-dmso)


##AP1-associated genes (Figure 3E)

source('../get_nearest_TE.R') ##hijacking this function to use with a peak file

ap1 = get_nearest_TE(rna.av, '../../Peaks/hTSC_CutTag_cJun.relaxed.bed') %>%
	mutate(is.near = distance<5000 & distance>=0) %>%
	filter(dmso>7 | sp>7)

quartz(w=1.5, h=3.5)
ggplot(ap1, aes(x=is.near, y=fc)) +
	geom_boxplot(outlier.shape=NA, fill='lightblue') +
	ylim(-2,3) +
	theme_classic() +
	xlab('') +
	ylab('log2 fold change')

wilcox.test(ap1$fc ~ ap1$is.near)


##TE-associated genes (Figure 3F)

te = get_nearest_TE(rna.av, '../../Annotations/H3K27ac_TEs.bed') %>%
	filter(distance<100000 & distance>=0) %>%
	filter(dmso>7 | sp>7)

quartz(w=5, h=3.5)
filter(te, !(te.name %in% c('MER11D','MER39B','MER41A','MER41C','MER61D'))) %>%
ggplot(aes(x=te.name, y=fc)) +
	geom_boxplot(outlier.shape=NA, fill='lightblue') +
	geom_hline(yintercept=0, linetype='dashed') +
	ylim(-2,3) +
	theme_classic() +
	xlab('') +
	ylab('log2 fold change')


wt = group_by(te, te.name) %>% summarise(p = wilcox.test(fc)$p.value)
p.adjust(wt$p, method='BH')
