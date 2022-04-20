
#Okae data

rna = read.delim('Okae_RNA_hTSC.txt.gz',as.is=T)


#Match to Ensemble IDs and hg38 coordinates via Biomart

bio = read.delim('Ensembl_hg38_biomart.txt.gz',as.is=T)

id.match = match(rna$Gene.symbol, bio$Gene.name)   #match primary gene names
sec.match = match(rna$Gene.symbol, bio$Gene.Synonym)   #match gene synonyms
id.match[is.na(id.match)] = sec.match[is.na(id.match)]   #merge

gene = bio$Gene.name[id.match]
id = bio$Gene.stable.ID[id.match]
chr = paste('chr',bio$Chromosome.scaffold.name[id.match],sep='')
start = bio$Gene.start..bp.[id.match]
end = bio$Gene.end..bp.[id.match]
strand = bio$Strand[id.match]
strand[strand==1] = '+'
strand[strand==-1] = '-'


##Make hg38 table with IDs and official gene symbols

df = data.frame(gene,id,chr,start,end,strand,rna[,3:ncol(rna)],stringsAsFactors=FALSE)

chr.list = c(paste('chr',c(1:22),sep=''),'chrX','chrY')
df.main = df[df$chr %in% chr.list,]   #only keep genes from assembled chromosomes


##Deal with duplicated genes
#selects the entry wiht highest mean expression across all samples

df.sort = df.main[order(df.main$gene),]

mean.expr = rowSums(df.sort[7:ncol(df.sort)])

dup = logical(nrow(df.sort))
for (i in 2:nrow(df.sort)) {
	max.id = i-1
	if (df.sort$gene[i]==df.sort$gene[i-1]) {
		if (mean.expr[i]>mean.expr[max.id]) {
			max.id = i
			dup[i-1] = TRUE
		} else {
			dup[i] = TRUE
		}
	}
}

dedup = df.sort[!dup,]


##Write output

write.table(dedup,'RNA_hTSC_hg38.txt',sep='\t',quote=F,row.names=F)

