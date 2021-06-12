#!/bin/bash 
###################################################
#
#	Bash script to generate heatmaps and annotation
#	tables form Revigo treemaps
#Author: Arturo Vera
###################################################

#activate R
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate R_env

script='/home/avera/scripts/Revigo/Aftermsystems'
for i in *.BP.R;
	do
	a=$(echo $i|sed 's/.R//g');
	b=$(echo $i|sed 's/.GOseq.enriched.revigo.BP.R//g');
	perl $script/mod.revigoR.iso.pl \
	$i \
	$a.collapsed.tab \
	$script/plantilla.mod.isoforms.txt \
	$a.heatmapB.pdf \
	$b.Cluster.Annotation.DEG.tab > $a.revigoheatmap.R;
	Rscript $a.revigoheatmap.R
done
