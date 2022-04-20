library(gplots)


##read data

ortho.files = list.files('./orthologues')

ortho = list()
for (i in 1:length(ortho.files)) {
	ortho[[i]] = read.delim(paste('./orthologues', ortho.files[i], sep='/'), as.is=T, header=F)
}
names(ortho) = gsub('_orthologues.bed','',ortho.files)


##count elements per family

families = scan('te_families.txt', character())

count = numeric()
hg38.count = numeric()
for (fam in paste(families,'_',sep='')) {
	count = rbind(count, unlist(lapply(ortho, function(x) length(grep(fam, x$V4)))))
	hg38.count = c(hg38.count, nrow(read.delim(paste('../Annotations/',fam,'hg38.bed',sep=''),as.is=T)))
}
rownames(count) = families
count = cbind(count, hg38.count)
colnames(count)[ncol(count)] = 'hg38'


##get percentage

perc = count/matrix(rep(hg38.count, ncol(count)), ncol=ncol(count)) * 100


##re-order and transpose

sp.order = c('hg38','panTro6','gorGor6','ponAbe3','nomLeu3','rheMac10','calJac4','tarSyr2','micMur2')
te.order = c("MER31A","MER39","LTR23","MER39B","MER21A","MER61E",
			"MER61D","MER41C","MER41B","MER41A","LTR8B","LTR8",
			"LTR10A","LTR10F","LTR7C","LTR3A","MER11D","LTR2B")

perc2 = t(perc[match(te.order, rownames(perc)), match(sp.order, colnames(perc))])


##heatmap (Figure 4B)

svg('~/Downloads/orthology.svg', width=6, height=4.5)
heatmap.2(perc2,
		scale='none',
		Rowv=F, Colv=F, dendrogram='none',
		trace='none',
		col=colorRampPalette(c('white','red')),
		colsep=1:ncol(perc2), rowsep=1:nrow(perc2), sepwidth=c(0.02, 0.02),
		density.info='none',
		breaks=seq(0,100,2))
dev.off()
