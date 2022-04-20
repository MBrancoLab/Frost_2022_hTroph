##Download Repeatmasker annotation and convert to bed files


genome = 'hg38'

##download annotations

system(paste('wget http://www.repeatmasker.org/genomes/',
	genome, '/RepeatMasker-rm405-db20140131/',
	genome, '.fa.out.gz', sep=''))


##read data

rmask = read.delim(paste(genome,'.fa.out.gz',sep=''),
	sep='', skip=3, as.is=T, header=F)


##make bed file

strand = rmask$V9
strand[strand=='C'] = '-'

bed = data.frame(chr=rmask$V5, start=rmask$V6, end=rmask$V7,
	name=rmask$V10, score=rmask$V1, strand)


##write full table without simple repeats (to use in permutation analysis)

bed.f = bed[rmask$V11!='Simple_repeat' & rmask$V11!='Low_complexity',]

write.table(bed.f, paste('Repeatmasker_', genome, '_4.0.5_noSimple.bed', sep=''),
	sep='\t', quote=F, row.names=F, col.names=F)


##get H3K27ac-associated elements from families of interest

sel = scan('../Active_families/te_families.txt', character())
bed.s = bed[bed$name %in% sel,]

write.table(bed.s, 'sel_families.temp', sep='\t', quote=F, row.names=F, col.name=F)
system(paste('~/Documents/bedtools2/bin/intersectBed -u -f 0.5 -a sel_families.temp -b ../Peaks/hTSC_CutTag_H3K27ac.relaxed.bed > ../Annotations/H3K27ac_TEs.bed', sep=''))
unlink('sel_families.temp')

