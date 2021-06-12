#!/bin/bash
##############################################################################
#	Goseq.sh Script obtain GO enrriched terms from Differential expressed
#	isoforms.
#	Author: Arturo Vera
#
#	Dependecies: Trinity utility script: run_GOseq.pl
#			Background.isoforms
#			isoform.length.txt
#			GO.transcripts.ancestry.tab
#			Union cluster tables from h.clust.R script
#			
#################################################################################

#generating Factor
scripts='/home/avera/scripts/GO'
ln -s $scripts/Background.isoforms .
ln -s $scripts/isoform.length.txt .
ln -s $scripts/GO.transcripts.ancestry.tab .

echo "Generating factors"
for i in Union*Cl*.tab;
	do
	a=$(echo $i|cut -d . -f 1,2);
	cat $i|perl -e '($file,$name)=@ARGV;open(FH, $file);while(<FH>){chomp; if($_ =~ /^TR/){@col=split(/\t/); print "$name\t$col[0]\n";}}' $i $a > $i.factor.txt ;

done

#Running GOSeq
echo "Running GOseq"
for i in *factor.txt;
	do
	/home/avera/bin/trinityrnaseq-v2.12.0/Analysis/DifferentialExpression/run_GOseq.pl --factor_labeling $i --GO_assignments GO.transcripts.ancestry.tab --lengths isoform.length.txt --background Background.isoforms

done
