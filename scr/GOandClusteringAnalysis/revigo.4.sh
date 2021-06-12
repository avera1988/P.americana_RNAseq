#!/usr/bin/bash
##########################################################################################
#	Script for Collaps Go enriched terms by REvigo
#	Author Arturo Vera
#	All dependencies are in ~/scripts/Revigo/Aftermsystems
#
##########################################################################################
if [ $# -eq 0 ]; then
	echo "usage: $0 enriched Sematical_Clustering_cutoff_REVIGO(0.9,0.7,0.5 or 0.4) Num_clusters";
	exit 1;
	fi

file=$1
revigoSC=$2
START=1;
END=$3
scripts='/home/avera/scripts/Revigo/Aftermsystems'

#Modifiying enriched go files For Revigo
echo "Preparing files"

for i in *$file ;
	do
	cat $i |awk '{if (NR!=1) {print}}' |cut -f 1,2 > $i.revigo;
done

echo "REVIGO"


#retrieving REVIGO it takes some time
for i in *.revigo ;
	do
	python $scripts/scrape_revigo.mod.html5.py $i $i $revigoSC;
done

#Run the R scripts for Treemaps

echo "Generating Treemaps"

#activate R
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate R_env

for i in *.R ; 
	do 
	name=$(echo $i|cut -d . -f 1,2,3,4,5)
	a=$(echo $i|cut -d . -f 6);
	#b=$(echo $i|cut -d . -f 1,2);
	Rscript $name.$a.R ; 
	mv revigo_treemap.pdf $name.$a.revigo_treemap.pdf;
done


echo "Collapsing and correlated groups"

#Collaps all redundant GO parents classified by REVIGO and correlate Trinity ID for Boxplots for each cluster
for ((i = $START ; i <= $END ; i++))
	do 
	echo "correlated";
	#cat REVIGO.C${i}.BP.csv| sed 's/"//g'  > REVIGO.C${i}.BP.csv.mod;
	Rscript $scripts/revigo.correlated.3.R Union.${i}.GOseq.enriched.revigo.BP.csv Union.${i}.GOseq.enriched; 
	#Rscript $scripts/revigo.correlated.2.R Union.${i}.GOseq.enriched.revigo.CC.csv Union.${i}.GOseq.enriched;
	#Rscript $scripts/revigo.correlated.2.R Union.${i}.GOseq.enriched.revigo.MF.csv Union.${i}.GOseq.enriched;
	echo "collapsing";
	cat Union.${i}.GOseq.enriched.revigo.BP.csv.correlated.tab | sed 's/"//g' > Union.${i}.GOseq.enriched.revigo.BP.csv.correlated.tab.mod;
	#cat Union.${i}.GOseq.enriched.revigo.CC.csv.correlated.tab | sed 's/"//g' > Union.${i}.GOseq.enriched.revigo.CC.csv.correlated.tab.mod;
	#cat Union.${i}.GOseq.enriched.revigo.MF.csv.correlated.tab | sed 's/"//g' > Union.${i}.GOseq.enriched.revigo.MF.csv.correlated.tab.mod;
	perl $scripts/collapsing.revigo.v.2.pl Union.${i}.GOseq.enriched.revigo.BP.csv.correlated.tab.mod > Union.${i}.GOseq.enriched.revigo.BP.collapsed.tab;
	#perl $scripts/collapsing.revigo.pl Union.${i}.GOseq.enriched.revigo.CC.csv.correlated.tab.mod > Union.${i}.GOseq.enriched.revigo.CC.collapsed.tab;
	#perl $scripts/collapsing.revigo.pl Union.${i}.GOseq.enriched.revigo.MF.csv.correlated.tab.mod > Union.${i}.GOseq.enriched.revigo.MF.collapsed.tab;
done

#Moving into a specific directories

mkdir BP_dir
#mkdir MF_dir
#mkdir CC_dir

mv *.BP.* BP_dir
#mv *.MF.* MF_dir
#mv *.CC.* CC_dir

