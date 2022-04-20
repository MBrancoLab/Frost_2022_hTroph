ucscDir=~/Documents/UCSC_tools

bed=$(sed s/$/_hg38.bed/ te_families.txt | tr '\n' ' ')

cd ../Annotations
cat $bed > ../TE_families/tsc_tes.bed
cd ../TE_families


for sp in panTro6 gorGor6 ponAbe3 nomLeu3 rheMac10 calJac4 tarSyr2 micMur2
do
	upSp="$(tr '[:lower:]' '[:upper:]' <<< ${sp:0:1})${sp:1}"
	${ucscDir}/liftOver.dms -minMatch=0.1 tsc_tes.bed ${ucscDir}/hg38To${upSp}.over.chain temp.bed failed.bed
	${ucscDir}/liftOver.dms -minMatch=0.1 temp.bed ${ucscDir}/${sp}ToHg38.over.chain orthologues/${sp}_orthologues.bed failed.bed
done
rm temp.bed failed.bed tsc_tes.bed
