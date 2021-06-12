#!/usr/bin/perl
######################################################################################################################################################	
#		This script merge the DE analysis Tables form DESeq wiht the Annotation Profile form each expressed gene and transcript
#	Author: Arturo Vera
# 	October 2018
#	
#	It needs The annotation of DE genes from Annotation.DE.sh, LFC and p-values from DESeq analsysis obtained from the same shell script
#	
#
########################################################################################################################################################

use strict;
use warnings;
#Defining variables
my($annot, $DE, $final, $depval, $id, $lfc, $pval, %filas,$Gen_ID,$Isoform_ID,$Uniprot_ID,$Uniprot_annot,$PFAM_ID,$PFAM_annot,$BLASTX_Invertebrate_ID,$BLASTX_Invertebrate_annot,$AMP, $id2,$pv,$anot2, %key, $gen_id,$isoform_id,$lc,$Uniprot_id,$uniprot_annot,$pfam_id,$pfam_annot,$blastx_invertebrate_id,$blastx_invertebrate_annot,$amp,$total,$gen_Id,$Isoform_id,$LC,$Pvalue,$Uniprot_Id,$uniprot_Annot,$Pfam_id,$pfam_Annot,$Blastx_invertebrate_id,$Blastx_invertebrate_annot,$Amp,$ko,$kegg,$cazy,$cazyan,$Ko,$Kegg,$Cazy,$Cazyan,$KO,$KEGG,$CAzy,$CAzyan);
@ARGV == 3 or die "Usage: perl $0 DE_Table_annotated\tDE_Table\tfinal_file_nale\n";
($annot, $DE, $final)=@ARGV;
open (DE, $DE);
open (ANNOT, $annot);
open DEPVAL, ">tmp.DE.pval.tab";
open TMP_ANOT, ">tmp.DE.annot.tab";
$depval='tmp.DE.pval.tab';
$anot2='tmp.DE.annot.tab';
open (ANOT2, $anot2 );
open (FH1, $depval);
open FINAL, ">$final.tab";
$total=$final.'.tab';
open (FH2, $total);
open UP, ">$final.up.tab";
open DOWN, ">$final.down.tab";
print "Staring to merge DE and annotation\n";
#Reading DE table
while(<DE>){
	chomp;
	if($_ =~ /^T/){
	($id, $lfc, $pval)=split(/\t/);
	print DEPVAL "$id\t$pval\n";
	$filas{$id}=$lfc;
	}
}
close DE;
close DEPVAL;
#Reading Trinotate annotation table
while(<ANNOT>){
	chomp;	
	if($_ =~ /^T/){	
	($Gen_ID,$Isoform_ID,$Uniprot_ID,$Uniprot_annot,$PFAM_ID,$PFAM_annot,$BLASTX_Invertebrate_ID,$BLASTX_Invertebrate_annot,$ko,$kegg,$cazy,$cazyan,$AMP)=split(/\t/);
	if(exists $filas{$Isoform_ID}){
	print TMP_ANOT "$Gen_ID\t$Isoform_ID\t$filas{$Isoform_ID}\t$Uniprot_ID\t$Uniprot_annot\t$PFAM_ID\t$PFAM_annot\t$BLASTX_Invertebrate_ID\t$BLASTX_Invertebrate_annot\t$ko\t$kegg\t$cazy\t$cazyan\t$AMP\n"; 
	}
 }
}
close(ANNOT);
close(TMP_ANOT);
#reading p-values
while(<FH1>){
	chomp;
		($id2, $pv)=split(/\t/);
		#print "$id2\n";
		$key{$id2}=$pv;
}
#Merging LFC, p-val and Annotation
print FINAL "Gene_id\tIsoform_id\tLFC\tp_value\tUniref_BX_ID\tUniref_Annot\tPFAM_ID\tPFAM_Annot\tInvertebrateID\tInvertebrate_Annotation\tKO\tKEGG_annot\tCAZyME\tCAZY_annotation\tAMP\n";
while(<ANOT2>){
		chomp;
		($gen_id,$isoform_id,$lc,$Uniprot_id,$uniprot_annot,$pfam_id,$pfam_annot,$blastx_invertebrate_id,$blastx_invertebrate_annot,$Ko,$Kegg,$Cazy,$Cazyan,$amp)=split(/\t/);
		if(exists $key{$isoform_id}){
			print FINAL "$gen_id\t$isoform_id\t$lc\t$key{$isoform_id}\t$Uniprot_id\t$uniprot_annot\t$pfam_id\t$pfam_annot\t$blastx_invertebrate_id\t$blastx_invertebrate_annot\t$Ko\t$Kegg\t$Cazy\t$Cazyan\t$amp\n";
		}
}
close(ANOT2);
close(FH1);
close (FINAL);
print "Generating Up and Down files\n";
#obtaining up and down
print UP "Gene_id\tIsoform_id\tLFC\tp_value\tUniref_BX_ID\tUniref_Annot\tPFAM_ID\tPFAM_Annot\tInvertebrateID\tInvertebrate_Annotation\tKO\tKEGG_annot\tCAZyME\tCAZY_annotation\tAMP\n";
print DOWN "Gene_id\tIsoform_id\tLFC\tp_value\tUniref_BX_ID\tUniref_Annot\tPFAM_ID\tPFAM_annot\tInvertebrateID\tInvertebrate_Annotation\tKO\tKEGG_annot\tCAZyME\tCAZY_annotation\tAMP\n";
while(<FH2>){
	chomp;
	if($_ =~ /^T/){
			($gen_Id,$Isoform_id,$LC,$Pvalue,$Uniprot_Id,$uniprot_Annot,$Pfam_id,$pfam_Annot,$Blastx_invertebrate_id,$Blastx_invertebrate_annot,$KO,$KEGG,$CAzy,$CAzyan,$Amp)=split(/\t/);
			#print "perro\n";
			if($LC < 0){
			print DOWN "$gen_Id\t$Isoform_id\t$LC\t$Pvalue\t$Uniprot_Id\t$uniprot_Annot\t$Pfam_id\t$pfam_Annot\t$Blastx_invertebrate_id\t$Blastx_invertebrate_annot\t$KO\t$KEGG\t$CAzy\t$CAzyan\t$Amp\n";	
	}
	if($LC >0){
		print UP "$gen_Id\t$Isoform_id\t$LC\t$Pvalue\t$Uniprot_Id\t$uniprot_Annot\t$Pfam_id\t$pfam_Annot\t$Blastx_invertebrate_id\t$Blastx_invertebrate_annot\t$KO\t$KEGG\t$CAzy\t$CAzyan\t$Amp\n";
	}
 }
}
close(UP);
close(DOWN);
close(FH2);
print "I have finished\n";
