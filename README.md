# Frost_2022_hTroph

Scripts and files used in:

Frost JM, Amante SM, Okae H, Jones EM, Ashley B, Lewis RM, Cleal JK, Caley MP, Arima T, Maffucci T, Branco MR (2022) "Regulation of human trophoblast gene expression by endogenous retroviruses" *bioRxiv*


## Primary processing

*ChIP-seq*

Data were aligned using the bowtie2.sh script, which also produces bigwig tracks. Peaks (in 'Peaks' folder) were detected using the macs2.sh script.

*CUT&Tag/CUT&RUN*

Data were aligned using the bowtie2.sh script, which also produces bigwig tracks. Peaks (in 'Peaks' folder) were detected using seacr_prep.sh followed by seacr.sh.

*RNA-seq*

Data were aligned using the hisat.sh script.


## Peak enrichment

*Permutation test*

A permutation test was used to identify TE families enriched for histone modifications and transcription factors. The main script is permutation.sh, which runs the permutation_hg38.R script for a list of peak files from either MACS2 or SEACR. See permutation_hg38.R for more details about other files used. The output files contain details of the number of overlapping peaks for each TE family, with associated stats.

*H3K27ac-enriched families*

The overlap files from the permutation test are used by the k27ac_enrichment.R script to identify H3K27ac-enriched TE families in hTSCs or hESCs, and produce the enrichment scatter plots in Figure 1B.

*ENCODE DHS*

The same permutation test as above is run by permutation_encode.sh (which wraps around permutation_encode.R) on a list of pre-processed bed files from ENCODE with peaks from DNase-seq data for different tissues. The concatenated enrichment values (DHS_enrichments.txt) are used by DHS_enrichment.R to produce the plots in Figure 1C and Supplementary Figure 1B. This also defines criteria to narrow down the list of TE families of interest.


## Annotations

The RMasker_outToBed.R script downloads the Repeatmasker file used here and can create bed annotations for any TE family of interest (several are stored within this directory). It can also produce an annotation of all H3K27ac-overlapped TEs from families of interest, using as input the final lists of candidate families from the 'Active_families' folder (te_families.txt for hTSC-active families; esc_families.txt for hESC-active ones) and the respective H3K27ac peaks from the 'Peaks' folder.


## Active families

*H3K27ac heatmaps*

The H3K27ac heatmaps in Figure 1C and Supplementary Figure 1A were produced with the heatmaps.sh script.

*Histone modifications*

The histone_mods.R script takes the list of hTSC-active TE families and peak files for different histone modifications to produce the circos plots in Figure 1E.

*EVT activity*

The evt_activity.R script first compares the H3K27ac enrichment values between hTSCs and EVT for selected TE families (producing Supplementary Figure 1D). It then intersects EVT-enriched families with several histone modification peaks to classify EVT-active TEs based on their hTSC chromatin state (producing Supplementary Figure 1F).

*TE orthology*

The liftover.sh script takes the list of hTSC-active TE families and runs reciprocal liftOver between the hg38 assembly and that of selected non-human primates (liftOver chain files available via UCSC). Lists of orthologous TEs from this script are in the 'orthologues' folder. The orthology.R script then takes this output and produces the heatmap in Figure 4B.


## Transcription factors

*Motif analysis*

The meme.sh script uses the MEME suite to first identify enriched motifs for each TE family of interest using the AME tool (output in 'AME' folder). It then takes the list of enriched motifs and runs FIMO to pinpoint these motifs in individual TE copies (output in 'FIMO' folder). A list of motifs of interest (selected_motifs.txt) was manually compiled based on outputs from select_motifs.R and literatures searches. The FIMO_summary.R script then this list of motifs and the FIMO output to produce the plots in Figure 2A and Supplementary Figure 2A. The alignments in Figure 2D were produced by the family_comparisons.R script.

*Binding enrichment*

The TF_enrichments.R takes the permutation test results from transcription factor profiling data and the lists of hTSC- and hESC-active families to produce the scatter plots in Figure 2C and Supplementary Figure 2B. Binding profiles for JunD and c-Jun across families of interest were produced by ap1_profiles.sh (output in 'AP1_profiles' folder) and ap1_profiles.R (Figure 2E).


## RNA-seq

*Data sources*

1. Our raw data from primary cytotrophoblast published in Ashley et al (PMID:35256050). These were used for identifying genes with TE-derived promoters.
2. Processed data from Okae et al (PMID:29249463) for primary and hTSC-derived trophoblast cells. Data were downloaded (Okae_RNA_hTSC.txt.gz) and coupled to gene names using tidy_Okae_RNAseq.R and biomart data (Ensembl_hg38_biomart.txt.gz), producing RNA_hTSC_hg38.txt.gz.
3. Processed data from Dong ett al (PMID:32048992) for hESC-to-hTSC transdifferentiation (Dong_RNA_transDiff.txt.gz).
4. Our data from hTSCs treated with the JNK inhibitor SP600125. Raw read counts were produced using Seqmonk (AP1_read_counts.txt file in JNK_inhibition folder).

*TE-derived promoters*

The stringtie.sh script uses Stringtie to perform transcriptome assembly, then takes the transcriptional start sites of multi-exonic transcripts (using get_TSS.R) and intersects it with the Repeatmasker annotation (output in te_promoters folder). The TE_promoters.R script then takes these transcripts and selects those associated with hTSC-active families, producing TE_promoters.txt.


*Expression in trophoblast*

The RNAseq_analysis.R script first renormalises the Okae data, then finds the active TEs nearest to each gene promoter (using get_nearest_TE.R) and plots their expression in trophoblast relative to stroma as a function of the TE-gene distance (Figure 3A and Supplementary Figure 3B). The script also plots the relative expression of genes close to active TEs in hTSCs or hESCs using the Dong dataset (Figure 3B). The heatmap in Figure 3C was produced by focusing on variably expressed genes (Okae data). Finally, this script breaks down the relative trophoblast expression (Okae data) by family (Supplementary Figure 3C).

*JNK inhibition*

The AP1_DEseq.R script performs differential expression analysis on the JNK inhibition data, and produces normalised gene expression values (AP1_RNA_vsd.txt). It then performs gene ontology analysis uwing topGO (Supplementary Figure 4A) and plots the expression of selected genes (Supplementary Figure 4B). The fold chance of JUN-associated gene expression is plotted based on gene promoter proximity to JUN peaks (Figure 3E). Finally, the expression fold change for genes close to hTSC-active TEs is plotted, split by family (Figure 3F).


## Comparative analysis

*Data pre-processing*

All RNA-seq data were mapped using the hisat.sh script in the 'Primary_processing' folder and RPKM values extracted using Seqmonk. Human data (hTSC_RNA_RPKM.txt.gz) is from Okae et al (PMID:29249463), macaque data (macTSC_RNA_RPKM.txt.gz) from Schmidt et al (PMID:33154556), and mouse data (mTSC_RNA_RPKM.txt.gz) from Cambuli et al (PMID:25423963). Mouse H3K27ac peaks (from CUT&Tag data processed as described above) were intersected with Repeatmasker annotations for RLTR13D5 and RLTR13B elements, producing mTSC_H3K27ac_TEs.bed.

*Orthologous gene expression*

The comparative_analysis.R script first parses the Seqmonk-generated RPKM tables, and adds information on the nearest active TE for human and mouse. It then merges the human data with either mouse or macaque data, based on 1-to-1 gene orthologues (orthologous_genes.txt.gz, from biomart), and performs a quantile-based normalisation. In the human-macaque comparison it also determines whether the nearby active TE is orthologous between the two species. It finally produces boxplots for the difference in gene expression between species based on the classification of the nearest active TE (Figures 4C and 4D).


## Assays

*RT-qPCR*

The RT-qPCR.R script imports GraphPad Prism files containing normalised data and produces the expression plots in Figures 5 and 6, and Supplementary Figures 5 and 6.

*FACS histograms*

The FACS.R script imports FCS-formatted files from FACS analysis and produces the histograms in Supplementary Figure 6B.

*Other data*

The other_data.R script imports GraphPad Prism files and produces: a) the growth curves in Figure 6C, b) the HLA-G+ FACS counts in Figure 6D), c) the sENG ELISA results in Figure 6E.



