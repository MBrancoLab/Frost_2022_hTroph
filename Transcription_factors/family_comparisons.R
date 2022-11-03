
library(msa)

clean.cons = function(consensus) {
	nodash = gsub('-','',consensus)
	return(gsub('\\.','N',nodash))
}



##LTR10

ltr10 = readDNAStringSet('fasta/LTR10_dfam.fa')
ltr10.mus = msa(ltr10, method='Muscle', type='dna')

ltr10a = readDNAStringSet('fasta/TSC_LTR10A.fa')
ltr10a.mus = msa(ltr10a, method='Muscle', type='dna')
ltr10a.cons = clean.cons(msaConsensusSequence(ltr10a.mus, type='upperlower', thresh=c(50,25)))

ltr10f = readDNAStringSet('fasta/TSC_LTR10F.fa')
ltr10f.mus = msa(ltr10f, method='Muscle', type='dna')
ltr10f.cons = clean.cons(msaConsensusSequence(ltr10f.mus, type='upperlower', thresh=c(50,25)))

ltr10.merge = c(ltr10[c(1,2,8,9)], ltr10a.cons,ltr10f.cons)
names(ltr10.merge)[5:6] = c('ts_ltr10a','ts_ltr10f')

print(msa(ltr10.merge,method='Muscle',type='dna'), show='complete')


##LTR8

ltr8.dfam = readDNAStringSet('fasta/LTR8_dfam.fa')

ltr8 = readDNAStringSet('fasta/TSC_LTR8.fa')
ltr8.mus = msa(ltr8, method='Muscle', type='dna')
ltr8.cons = clean.cons(msaConsensusSequence(ltr8.mus, type='upperlower', thresh=c(50,25)))

ltr8b = readDNAStringSet('fasta/TSC_LTR8B.fa')
ltr8b.mus = msa(ltr8b, method='Muscle', type='dna')
ltr8b.cons = clean.cons(msaConsensusSequence(ltr8b.mus, type='upperlower', thresh=c(50,25)))

ltr8.merge = c(ltr8.dfam, ltr8.cons, ltr8b.cons)
names(ltr8.merge)[4:5] = c('ts_ltr8','ts_ltr8b')

print(msa(ltr8.merge,method='Muscle',type='dna'), show='complete')


