#!/usr/bin/perl
use strict;
my(@line,$geneid,%seen);
print "Gene_ID	transcript_id\tUniref_BX_ID\tUniref_Annot\tPFAM_ID\tPFAM_Annot\tInvertebrateID\tInvertebrate_Annotation\tKO_number\tKEEG_annot\tCAZyME\tCAZY_annotation\tAMP\n";
while(<>){
	chomp; 
	@line=split(/\t/); 
	$geneid=$line[0]; 
	print "$geneid\t$line[1]\t$line[2]\t$line[3]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$line[10]\t$line[11]\t$line[12]\t$line[13]\n"  if ! $seen{$line[1]}++;
}
