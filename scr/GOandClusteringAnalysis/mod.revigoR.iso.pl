#!/usr/bin/perl
#######################################################################
#	This script parses the R script for treemap from Revigo
#	and generates a heatmap and tables with annotation of each
#	transcript.
#	Dependencies:
#		- RevigoTreemap.Rfile; 
#		- Collapsed table from Revigo
#		- plantilla.isoforms.txt for R script generation
#		- Annotation and DEG values in each cluster
#	Author: Arturo Vera
#######################################################################
use strict;

@ARGV == 5 || die "Usage: $0 RevigoTreemap.Rfile collapsed.tab plantilla.mod.isoforms.txt heatmap.pdf .Annotation.DEG.tab\n";

my($rfile,$datos,$plant,$heatmapB,$annot)=@ARGV;
my(@rfile,$lines);
open(RFILE,$rfile);
open(PLAN,$plant);


while(<RFILE>){
	chomp;
	if($_ =~ /^treemap/){ 
	print "tm <- $_\n";
	}elsif($_ == /^pdf/ & $_ == /^dev/){
		$lines=$_;
		print "$lines\n";
	#}elsif($_ ==/^treemap/){
	#	print "tmp <- $_";
	}
}
while(<PLAN>){
		if($_ =~ /^colapsed/){
			chomp;
			print "$_ read.delim\(\"$datos\",sep=\"\\t\",header=F\)\n";
		}elsif($_ =~ /^Annotation/){
			chomp;
			print "$_ read.delim\(\"$annot\",sep=\"\\t\",header=T\)\n";
	}else{
			print $_;
	}
}

print "ggsave\(phmapB, file=\"$heatmapB\",width = 12,height = 12\)\n";
#print "ggsave\(phmapT, file=\"$heatmapT\",width = 12,height = 12\)\n";
print "write.table(stuffExDEG.BAn,\"$annot.Revigo.tab\",sep=\"\\t\",row.names=F,quote=F)\n";

close(RFILE);
close(PLAN);
