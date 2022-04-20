library(flowCore)
library(tidyverse)


##data examples

ct = read.FCS('Specimen_001_CT +ve_001.fcs', transformation=FALSE)
ev = read.FCS('Specimen_001_EVT+ve_003.fcs', transformation=FALSE)
no.ab = read.FCS('Specimen_001_EVT-ve_004.fcs', transformation=FALSE)


##plot HLA-G histograms (Supplementary Figure 6B)

plot_hlag = function(fcs) {
	data = as_tibble(fcs@exprs)
	colnames(data) = c('Time','FCS.A','SSC.A','SSC.W','GFP','empty','BFP','HLA.G')
	
	ggplot(data, aes(x=log10(HLA.G+221))) +
		geom_histogram(aes(y=(..count..)/sum(..count..)*100), fill='lightblue', binwidth=0.02) +
		geom_vline(xintercept=2.8, linetype='dashed') +
		xlim(2, 4) +
		ylim(0, 11.5) +
		theme_classic() +
		xlab('HLA-G') +
		ylab('% Cells')
}

quartz(w=3, h=3)
plot_hlag(ct)
plot_hlag(ev)
plot_hlag(no.ab)

