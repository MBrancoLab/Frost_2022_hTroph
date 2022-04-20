#################################
#Annotate with nearest active TE#
#################################


get_nearest_TE = function(rna, te.file, family=NA) {

	##Make bed with gene starts	
	options(scipen=10)
	tss = rna$start
	tss[rna$strand=='-'] = rna$end[rna$strand=='-']
	tss.bed = data.frame(chr=rna$chr,start=tss,end=tss,stringsAsFactors=FALSE)
	
	##Sort data by gene starts
	tss.sort = tss.bed[order(tss.bed$chr,tss.bed$start),]
	rna.sort = rna[order(tss.bed$chr,tss.bed$start),]
	write.table(tss.sort,'tss_temp.bed',sep='\t',quote=F,row.names=FALSE,col.names=FALSE)
	
	##Sort TE file
	system(paste('~/Documents/bedtools2/bin/sortBed -i', te.file, '> te_sorted_temp.bed'))
	if (!is.na(family)) {
		file.rename('te_sorted_temp.bed', 'te_sorted_all_temp.bed')
		system(paste("grep $'\t'", family, "$'\t' te_sorted_all_temp.bed > te_sorted_temp.bed", sep=''))
		unlink('te_sorted_all_temp.bed')
	}

	##Bedtools closest	
	system('~/Documents/bedtools2/bin/closestBed -d -a tss_temp.bed -b te_sorted_temp.bed > out_temp.bed')
	closest = read.delim('out_temp.bed',as.is=T,header=F)
	unlink(c('tss_temp.bed','te_sorted_temp.bed','out_temp.bed'))

	##Couple data	
	rna.sort$te.name = closest$V7
	rna.sort$te.chr = closest$V4
	rna.sort$te.start = closest$V5
	rna.sort$te.end = closest$V6
	rna.sort$te.strand = closest$V9
	rna.sort$distance = closest$V10
	
	##output
	return(rna.sort)
}