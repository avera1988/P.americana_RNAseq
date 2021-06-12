#!/bin/bash 

script='/home/avera/scripts/Revigo'
for i in *.BP.R;
	do
	a=$(echo $i|sed 's/.R//g');
	b=$(echo $i|sed 's/.GOseq.enriched.revigo.BP.R//g');
	perl $script/mod.revigoR.iso.pl $i $a.collapsed.tab $script/plantilla.mod.isoforms.txt $a.heatmapB.pdf $b.Cluster.Annotation.DEG.tab\
	|sed 's/tmPlot/tm<-tmPlot/' > $a.revigoheatmap.R;
	Rscript $a.revigoheatmap.R
done
