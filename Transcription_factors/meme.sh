##Requires MEME suite
##Requires HOMER (bedtools can also be used for fasta extraction)

##AME; finds enriched motifs for all TE families
while read te; do

	#extract sequences from active copies
	grep ${te}$'\t' ../Annotations/H3K27ac_TEs.bed > $te-active.bed
	homerTools extract $te-active.bed ~/Documents/homer/data/genomes/hg38 -fa > $te-active.fa

	#run AME
	ame --control '--shuffle--' $te-active.fa ~/Documents/motif_databases/JASPAR/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt
	grep MA ame_out/ame.tsv | cut -f 3 >> motifs.txt

	#tidy up
	mv ame_out/ame.tsv AME/$te.tsv
	rm -r ame_out

done < ../Active_families/te_families.txt


##FIMO; gets TEs with motifs from all of the enriched above
sort motifs.txt | uniq > unique_motifs.txt
while read te; do

	#run FIMO for enriched motifs
	mopt=""
	while read motif; do
		mopt="$mopt --motif $motif"
	done < unique_motifs.txt
	fimo $mopt ~/Documents/motif_databases/JASPAR/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt $te-active.fa

	#tidy up
	mv fimo_out/fimo.tsv FIMO/$te.tsv
	rm -r fimo_out $te-active.bed $te-active.fa

done < ../Active_families/te_families.txt


##FIMO on inactive TEs and selected motifs
while read te; do

	#get inactive TEs
	intersectBed -v -a ../Annotations/${te}_hg38.bed -b ../Peaks/hTSC_CutTag_H3K27ac.relaxed.bed > $te-inactive.bed
	homerTools extract $te-inactive.bed ~/Documents/homer/data/genomes/hg38 -fa > $te-inactive.fa

	#run FIMO for selected motifs
	id=$(cut -f 1 selected_motifs.txt)
	mopt=""
	for motif in $id; do
		mopt="$mopt --motif $motif"
	done
	fimo $mopt ~/Documents/motif_databases/JASPAR/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt $te-inactive.fa

	#tidy up
	mv fimo_out/fimo.tsv FIMO_inactive/$te.tsv
	rm -r fimo_out $te-inactive.bed $te-inactive.fa

done < ../Active_families/te_families.txt