
library(tidyverse)

##functions

parse_file = function(file) {
	df = read.delim(file, row.names=1, skip=1, header=F)[,-1]
	df = df[,!is.na(df[2,])]
	tibble(bins=rep(as.numeric(df[1,]),3),
		rpm=c(as.numeric(df[2,]),as.numeric(df[3,]),as.numeric(df[4,])),
		tf=rep(c('cJun','JunD','IgG'),each=ncol(df)))
}

plot_trend = function(data) {
	ggplot(data, aes(x=bins, y=rpm, color=tf)) +
		theme_classic() +
		ylab('RPM') +
		labs(color='') +
		geom_line(size=1) +
		scale_colour_manual(values=c('lightblue','grey','darkblue')) +
		scale_x_continuous(name='', breaks=c(0,median(data$bins),max(data$bins)), labels=c('-3kb','centre','+3kb'))
}


##plot all families

pfiles = list.files('./AP1_profiles/')

for (pfile in pfiles) {
	profile = parse_file(paste('AP1_profiles', pfile, sep='/'))
	te = gsub('-profile.txt','',pfile)
	
	quartz(w=3.5,h=2.5)
	Sys.sleep(1)
	print(plot_trend(profile) + labs(title=te))
}


##plot selected (Figure 2E)

quartz(w=3.5,h=2.5)
plot_trend(parse_file('AP1_profiles/LTR10A-profile.txt'))
plot_trend(parse_file('AP1_profiles/LTR10F-profile.txt'))
plot_trend(parse_file('AP1_profiles/LTR10B-profile.txt')) + ylim(0,1.5)
plot_trend(parse_file('AP1_profiles/LTR8B-profile.txt')) + ylim(0,1.5)
plot_trend(parse_file('AP1_profiles/LTR8-profile.txt')) + ylim(0,1.5)
plot_trend(parse_file('AP1_profiles/LTR8A-profile.txt')) + ylim(0,1.5)