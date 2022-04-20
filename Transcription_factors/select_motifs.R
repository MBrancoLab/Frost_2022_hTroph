library(GenomicRanges)


##read Okae RNA-seq data

rna = read.delim('../RNA-seq/RNA_hTSC_hg38.txt.gz', as.is=T)


##read FIMO output

data = read.delim('FIMO/MER61E.tsv', as.is=T)


##find frequent motifs

m.freq = table(data$motif_alt_id)
frequent = names(m.freq[m.freq>0.5*length(unique(data$sequence_name))])


##find expressed motifs

genes = strsplit(gsub('\\([[:print:]]+\\)', '', names(m.freq)), split='::')
is.expr = unlist(lapply(genes, function(x) {
	expr = rna[rna$gene %in% toupper(x), grepl('^CT..', colnames(rna)) | grepl('^TS', colnames(rna))]
	if (any(apply(expr>1, 1, all))) {
		return(TRUE)
	} else {
		return(FALSE)
	}
}))
expressed = names(m.freq[is.expr])


##filter motifs

data2 = data[data$motif_alt_id %in% frequent & 
	data$motif_alt_id %in% expressed & 
	data$q.value<0.05,]


#####################################
# Cluster motifs to find redundancy #
#####################################


##make GenomicRanges object

gr = makeGRangesFromDataFrame(data2, seqnames.field='sequence_name')
gr$motif = data2$motif_alt_id


##initiate count matrix for overlaps between motifs

m.id = sort(unique(data2$motif_alt_id))
motifs = matrix(nrow=length(m.id), ncol=length(m.id))
motifs [,]= 0
rownames(motifs) = colnames(motifs) = m.id


##populate matrix by overlapping GenomicRanges subsets

for (m in m.id) {
	sub = gr[gr$motif==m]
	overlap = subsetByOverlaps(gr, sub, minoverlap=4)
	counts = table(overlap$motif)
	motifs[rownames(motifs)==m, match(names(counts),colnames(motifs))] = counts
}


##plot heatmap

heatmap(motifs, scale='none', col=hcl.colors(max(motifs)),
	cexRow=0.3, cexCol=0.3)

