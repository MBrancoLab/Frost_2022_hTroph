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
