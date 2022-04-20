
####################
#### Arguments  ####
####################

##bed file of regions of interest, e.g., peaks
args = commandArgs(trailingOnly=TRUE)
npf.file = args[1]

##RepeatMasker annotation (a subset with TE families from enrichment analysis only)
rm.file = 'rmask_subset.bed'

##file with chromosome sizes (from UCSC)
chrom.file = '~/genomes/GRCh38/hg38.chrom.sizes'

##optional: bed file of unmappable regions
unmap.file = '~/genomes/GRCh38/hg38_150PE_lowMap.bed.gz'

##path to bedtools
bedtools.path = ''



####################
#### Functions  ####
####################


##Make random control probes
##unmap.file is a bed file of unmappable regions of the genome
##To shuffle fragments anywhere in the genome, make unmap.file=''

make.random = function(bed,chrom.file,unmap.file='',seed=NA,bedtools.path='') {
	i.file = tempfile()
	write.table(bed,file=i.file,quote=F,sep="\t",col.names=F,row.names=F)
	
	out.file = tempfile()
	
	if (is.na(seed)) {
		seed = round(runif(1,0,1000000))
	}
	
	if (unmap.file=='') {
		command=paste(bedtools.path,"shuffleBed -i ", i.file," -g ",chrom.file," -seed ",seed," -noOverlapping > ",out.file,sep="")
	} else {
		command=paste(bedtools.path,"shuffleBed -i ",i.file," -g ",chrom.file," -seed ",seed," -excl ",unmap.file," -f 0 -noOverlapping > ",out.file,sep="")
	}	
	cat(command,"\n")
	try(system(command))
	res=read.table(out.file,header=F,as.is=T)
	
	unlink(c(i.file,out.file))
	return(res)
}


##intersectBED with Repeatmasker

intersectBED = function(bed,rm.file,opt.string="-u",bedtools.path='') {
	b.file = tempfile()
	write.table(bed,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
	out = tempfile()
	
	command = paste(bedtools.path,"intersectBed -a ",rm.file," -b ",b.file,' ',opt.string," > ",out,sep="")
	cat(command,"\n")
	try(system(command))
	res = read.table(out,header=F,as.is=T)

	unlink(c(b.file,out))
	return(res)
}



###############
#### Main  ####
###############


##read bed file

npf = read.delim(npf.file,header=F,as.is=T)
bed = npf[npf$V3-npf$V2<40000,1:3]  ##exclude enormous peaks


##get repeatmasker counts on selected families

rmask = read.delim(rm.file,header=F)
rm.count = tapply(rmask$V1,rmask$V4,length)


##intersect bed with repeatmasker

bed.rm = intersectBED(bed,rm.file,bedtools.path=bedtools.path)
real.overlap = tapply(bed.rm$V4,factor(bed.rm$V4,levels=levels(rmask$V4)),length)
real.overlap[is.na(real.overlap)] = 0


##shuffle bed and intersect

n=1000
rand.overlaps = matrix(nrow=length(real.overlap),ncol=n)
for (i in 1:n) {
	bed.random = make.random(bed,chrom.file=chrom.file,unmap.file=unmap.file,bedtools.path=bedtools.path)
	bed.rand.rm = intersectBED(bed.random,rm.file=rm.file,bedtools.path=bedtools.path)
	rand.overlaps[,i] = tapply(bed.rand.rm$V4,factor(bed.rand.rm$V4,levels=levels(rmask$V4)),length)
	rand.overlaps[is.na(rand.overlaps[,i]),i] = 0
}


##get enrichment

random.mean = rowMeans(rand.overlaps)
enrich = real.overlap/random.mean


##output results

out = data.frame(enrich)
write.table(out,paste(npf.file,'enrich.txt',sep='_'),sep='\t',quote=F,col.names=NA)

