library(tidyverse)


##get arguments

args = commandArgs(trailingOnly=TRUE)
gtf.file = args[1]
out.file = args[2]


##read gtf file

gtf = read_tsv(gtf.file, comment='#',
	col_names = c('chr','source','feature','start','end','score','strand','frame','attribute'))


##get transcript IDs

split.att = strsplit(gtf$attribute, split='[;\"]')

get.att = function(att) {
	unlist(lapply(split.att, function(x) {
			pos = grep(att, x)
			if (length(pos)>0) {
				return(x[pos + 1])
			} else {
				return('')
			}
		}))
}

gtf = add_column(gtf, transcript = get.att('transcript_id'),
					fpkm = get.att('FPKM'),
					gene = get.att('ref_gene_name'))


##filter multiexonic transcripts

multi.ex = gtf %>% filter(feature=='exon') %>% group_by(transcript) %>%
	summarise(n=length(transcript)) %>% filter(n>1)
	
multi.tr = gtf %>% filter(feature=='transcript', transcript %in% multi.ex$transcript) %>%
	mutate(tss = case_when(strand=='+' ~ start, strand=='-' ~ end))


##write bed file

tr.bed = multi.tr %>% select(chr, tss, gene, fpkm, strand) %>%
	mutate(end=tss+1, .after=tss)
write_tsv(tr.bed, out.file, col_names=FALSE)

