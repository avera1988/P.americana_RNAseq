#!/bin/bash
########################################################################################################
#	This script add annotations to DE genes list from DESeq2
#	It needs the annoations from Trinotate and other databases (see mod_trinotate.pl) 
#	It pars and generates DE tables with annotation Up and Down regulated genes and their fasta files
#	Author Arturo Vera
#	October 2018
#	Modification November 2018
#	Dependencies:
#	parseo_all.de.2.pl
#	merging_annotation.DE.5.pl
#	Trinotate.invertebrate.amp.kegg.cazy.tab
#	trinity.mod_header.fasta
########################################################################################################
if [ $# -eq 0 ]; then
        echo "usage: $0 <DE.Result.tab>";
        exit 1;
        fi


file=$1
scripts='/home/avera/scripts/Annotation_RNASeq'
ln -s $scripts/GO.ancestry.RData .
#Parsing DE tables
#Annotating
cat $file |awk '{if (NR!=1) {print}}'|cut -f 1 |fgrep -w -f - $scripts/Trinotate.invertebrate.amp.kegg.cazy.tab|perl $scripts/parseo_all.de.3.pl > $file.annot.trinotate.total.tab
cat $file.annot.trinotate.total.tab |cut -f 1|fgrep -w -v -f - $file > $file.no_annot.tab
#Getting LFC and P-values
cat $file |awk '{if(NF != 1){print}}'|awk '{print $1"\t"$3"\t"$NF}' > $file.mod


#Merging DE values and Annotation
perl $scripts/merging_annotation_DE.5.pl $file.annot.trinotate.total.tab $file.mod $file.annot.final

#sorting Up and Down
cat $file.annot.final.up.tab|(read -r; printf "%s\n" "$REPLY"; sort -r -nk3) > $file.annot.final.up.sorted.tab
cat $file.annot.final.down.tab|(read -r; printf "%s\n" "$REPLY"; sort -r -nk3) > $file.annot.final.down.sorted.tab

#Extracting Fasta from up and Down
cat $file.annot.final.up.tab|awk '{if (NR!=1) {print}}'|cut -f 2|fgrep -w -f - -A 1 $scripts/trinity.all.paried.alltreat.mapped.mod_header.fasta|grep -v '\-\-' > FASTA.$file.annot.final.up.fasta
cat $file.annot.final.down.tab |awk '{if (NR!=1) {print}}'|cut -f 2|fgrep -w -f - -A 1 $scripts/trinity.all.paried.alltreat.mapped.mod_header.fasta|grep -v '\-\-' > FASTA.$file.annot.final.down.fasta
#obtaing Stadistics
total_genes=$(cat $file |awk '{if (NR!=1) {print}}'|wc -l)
annot_genes=$(cat $file.annot.trinotate.total.tab|awk '{if (NR!=1) {print}}'|wc -l)
no_annot_genes=$(cat $file.no_annot.tab|awk '{if (NR!=1) {print}}'|wc -l)
up_genes=$(cat $file.annot.final.up.tab|awk '{if (NR!=1) {print}}'|wc -l)
down_genes=$(cat $file.annot.final.down.tab|awk '{if (NR!=1) {print}}'|wc -l)
echo "Number of total DE genes" $total_genes
echo "Number of Annotated DE genes" $annot_genes
echo "Number of Up-regulated genes" $up_genes
echo "Number of Down-regulated genes" $down_genes
echo "Number of non annotated DE genes" $no_annot_genes
#Adding GO
filemod=$(echo $file|sed 's/.tab//g')
Rscript $scripts/annot.go.R $file.annot.final.tab $filemod.annot.csv $file.annot.final.up.sorted.tab $filemod.up.csv $file.annot.final.down.sorted.tab $filemod.down.csv
unlink GO.ancestry.RData


