#########################################################################################################################################################
#DE_analysis using DEseq2 for multiple comparisons
#Author Arturo Vera
#date: Sep, 2018
#Modificated: Jan, 2019
#Description: This script runs DESeq2 for differential expression analysis quantification on Salmon transcripts counts.
#	it depends on DESeq2 (https://bioconductor.org/packages/release/bioc/html/DESeq.html) and tximport (http://bioconductor.org/packages/release/bioc/html/tximport.html) for count table preparation.
#	Useful notes : http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#input-data
#
#Example in R:
#	source("DESeq2_multigroup.R")
#	
#	Note: 
#	This scripts needs the "quant" directories from Salmon, also needs the transcripts to gen correlation table (cm_tx2gene.txt). It uses these two elements
#	to create the count matrix table. So be sure that you have these in the path you are running the script.
#	This will return:
#	PCA plots
#	Number of Counts
#	raw and normalized
#	Tables of DESeq2 results with out filtering 
#	Tables of differential expresed genes Total and above LFC 1.5
#	MA plots of each comparison
#	
########################################################################################################################################################
#Parallel with 8 cores, usually this script is running @ unity.asc.osu.edu whit 8 cpus and asking at least 20G RAM (qsub -I -X -l nodes=1:ppn=8 -l mem=20GB -l walltime=12:00:00)
library("BiocParallel")
register(MulticoreParam(16))
#Dependencies 

library(DESeq2)
library(tximport)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(EnhancedVolcano)	

#Reading files in directory
print (c("you are in:", getwd()))
print ("Geting directroy")
dir <- "./"
system("ls -l |grep quant|awk '{print $9}' > samples.id.txt")
run <- readLines("samples.id.txt")
files <- file.path(dir, run,"quant.sf")
names(files) <- run
print (c("All files exist:",all(file.exists(files))))
#importing genes and transcripts corelations
print("Starting files import")
tx2gene <- read_tsv("cm_tx2gene.txt")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
colnames(txi$counts)<- factor(sub("_transcripts_quant", "", colnames(txi$counts)))
#generating raw table
write.table(txi$counts, file="raw_Counts.tab", sep="\t", quote=F)
print ("Raw Counts matrix is on: raw_Counts.tab")
print ("Creating Sample Table")
sampleTable <- data.frame(condition=factor(sub("\\d", "",colnames(txi$counts))))
rownames(sampleTable) <- colnames(txi$counts)
colnames(sampleTable) <- "library"
sampleTable <- separate(sampleTable, library, into = c("treatment", "tissue"), sep = "_", remove = FALSE, extra = "drop")
sampleTable$tissue <- gsub("H", "hindgut", sampleTable$tissue)
sampleTable$tissue <- gsub("M", "midgut", sampleTable$tissue)
#creating DESeq objectt
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ treatment + tissue)
dds <- DESeq(dds, fitType = "local", parallel=TRUE)
#creating contrast: If you just want to compare treatment within each tissue, then the DESeq authors recommend to create a grouping variable by combining your two factors (check:https://gist.github.com/Perugolate/f0050430f08afac19ccedd827dc98d92)
sampleTable_i <- unite(sampleTable, treatment_tissue, c(treatment, tissue), remove = FALSE)
ddi <- DESeqDataSetFromTximport(txi, sampleTable_i, ~ treatment_tissue)
ddi <- DESeq(ddi, fitType = "local", parallel=TRUE)
conditionsName <- unique(as.vector(ddi$treatment_tissue))
#Obtaining fpkm
FPKM <- fpkm(ddi)
FpkmDF <- as.data.frame(FPKM)
save(FpkmDF, file="fpkm.RData")
#Obtainign Normalized counts DESEQ
table_counts_normalized <- counts(ddi, normalized=TRUE)
save(table_counts_normalized, file="DDS.Normalized.counts.RData")
print ("Starting Contrast comparison")
print ("Start Hindgut comparisons")
#Comparing contras for DE in hindgut (All odd numbers in ConditionsName vector)
##Start F_hindgut vs all others (GN and WT)

for (i in seq(3,5,2)){
	print (c(conditionsName[1],"vs",conditionsName[(i)]))
	#As we will use a p-adjust value cut-off <= 0.05 DESeq authors recomend an indepent filtering alpha=0.05 
	res05 <- results(ddi, contrast=c("treatment_tissue",conditionsName[1],conditionsName[(i)]),alpha=0.05)
	res05Ordered <- res05[order(res05$pvalue),]
	data05 <- na.omit(as.data.frame(res05Ordered))
	#writing all results 
	write.table(res05Ordered, paste0("DESEQ.",conditionsName[1],"_vs_",conditionsName[(i)],".tab"), sep="\t", quote=F)
	print(summary(res05))
	#Getting differential expressed genes below p-adjust value 0.05
	keepAllDEgenes <- (data05$padj<=0.05)
	genesDE<-data05[keepAllDEgenes,]
	print (c("Number of all DE genes:", dim(genesDE)[1]))
	write.table(genesDE, paste0("All.DE.genes.",conditionsName[1],"_vs_",conditionsName[(i)],".tab"), sep="\t", quote=F)
	#Obtainig genes with logFoldchage 1.5
	##keepAllDEgenesLFC2 <- (data05$padj<=0.05 & abs(data05$log2FoldChange) > 1.5)
	##genesDELFC2<-data05[keepAllDEgenesLFC2,]
	#new way for filtering using dplyr
	DF05 <- as.data.frame(res05)
	DF05$Gene_id <- row.names(res05)
	DF05 <- DF05[,c(7,1:6)]
	genesDELFC1.5 <- dplyr::filter(DF05, padj <= 0.05 & abs(log2FoldChange) >= 0.585)
	print (c("Number of DE genes with LFC_1.5:", dim(genesDELFC1.5)[1]))
	write.table(genesDELFC1.5, paste0("All.DE.genes.LFC1.5.",conditionsName[1],"_vs_",conditionsName[(i)],".tab"), sep="\t", quote=F,row.names=F)
	#Generating MA plots
	pdf(paste0("MA_plot.",conditionsName[1],"_vs_",conditionsName[(i)],".pdf"))	
	plotMA(results(ddi, contrast=c("treatment_tissue",conditionsName[1],conditionsName[(i)]),alpha = 0.05), main=paste0("DE Genes from:", conditionsName[1],"_vs_",conditionsName[(i)]))
	abline(h=c(-1.5,1.5),col="dodgerblue",lwd=2)	
	dev.off()
	#Generating Volcano plots with Enhanced Volcano (https://bioconductor.org/packages/devel/bioc/html/EnhancedVolcano.html)
	volcan <- EnhancedVolcano(res05,lab = "",x = "log2FoldChange",y = "pvalue", FCcutoff = 1.5,transcriptPointSize = 1.5,transcriptLabSize = 3.0,
        title = paste0("DE Genes from:", conditionsName[1],"_vs_",conditionsName[(i)]))
	ggsave(volcan, file=paste0("volcano.plot.",conditionsName[1],"_vs_",conditionsName[(i)],".png"))
}
##GN_hindgut vs GN_hindgut
print (c(conditionsName[3],"vs",conditionsName[5]))
#As we will use a p-adjust value cut-off <= 0.05 DESeq authors recomend an indepent filtering alpha=0.05 (https://www.pnas.org/content/107/21/9546.long)
res05 <- results(ddi, contrast=c("treatment_tissue",conditionsName[3],conditionsName[5]),alpha=0.05)
res05Ordered <- res05[order(res05$pvalue),]
data05 <- na.omit(as.data.frame(res05Ordered))
#writing all results
write.table(res05Ordered, paste0("DESEQ.",conditionsName[3],"_vs_",conditionsName[5],".tab"), sep="\t", quote=F)
print(summary(res05))
#Getting differential expressed genes below p-adjust value 0.05
keepAllDEgenes <- (data05$padj<=0.05)
genesDE<-data05[keepAllDEgenes,]
print (c("Number of all DE genes:", dim(genesDE)[1]))
write.table(genesDE, paste0("All.DE.genes.",conditionsName[3],"_vs_",conditionsName[5],".tab"), sep="\t", quote=F)
#Obtainig genes with logFoldchage 1.5
##keepAllDEgenesLFC2 <- (data05$padj<=0.05 & abs(data05$log2FoldChange) >= 1.5)
##genesDELFC2<-data05[keepAllDEgenesLFC2,]
#new way for filtering using dplyr
DF05 <- as.data.frame(res05)
DF05$Gene_id <- row.names(res05)
DF05 <- DF05[,c(7,1:6)]
genesDELFC1.5 <- dplyr::filter(DF05, padj <= 0.05 & abs(log2FoldChange) >= 0.585)
print (c("Number of DE genes with LFC_1.5:", dim(genesDELFC1.5)[1]))
write.table(genesDELFC1.5, paste0("All.DE.genes.LFC1.5.",conditionsName[3],"_vs_",conditionsName[5],".tab"), sep="\t", quote=F,row.names=F)
#Generating MA plots
pdf(paste0("MA_plot.",conditionsName[3],"_vs_",conditionsName[5],".pdf"))	
plotMA(results(ddi, contrast=c("treatment_tissue",conditionsName[3],conditionsName[5]),alpha=0.05), main=paste0("DE Genes from:", conditionsName[3],"_vs_",conditionsName[5]))
abline(h=c(-1.5,1.5),col="dodgerblue",lwd=2)	
dev.off()
volcan <- EnhancedVolcano(res05,lab = "",x = "log2FoldChange",y = "pvalue", FCcutoff = 1.5,transcriptPointSize = 1.5,transcriptLabSize = 3.0,
        title = paste0("DE Genes from:", conditionsName[1],"_vs_",conditionsName[(i)]))
ggsave(volcan, file=paste0("volcano.plot.",conditionsName[3],"_vs_",conditionsName[5],".png"))
print ("Finishing Hindgut comparisons")
print ("Starting Midgut comparisons")
#Comparing contrast for DE in midgut (All even numbers in conditionsName vector)
##Start FM vs GNM and WTM
for (i in seq(4,6,2)){
	#As we will use a p-adjust value cut-off <= 0.05 DESeq authors recomend an indepent filtering alpha=0.05 (https://www.pnas.org/content/107/21/9546.long)
	print (c(conditionsName[2],"vs",conditionsName[(i)]))
	res05 <- results(ddi, contrast=c("treatment_tissue",conditionsName[2],conditionsName[(i)]),alpha=0.05)
	res05Ordered <- res05[order(res05$pvalue),]
	data05 <- na.omit(as.data.frame(res05Ordered))
	#writing all results
	write.table(res05Ordered, paste0("DESEQ.",conditionsName[2],"_vs_",conditionsName[i],".tab"), sep="\t", quote=F)
	print(summary(res05))
	#Getting differential expressed genes below p-adjust value 0.05
	keepAllDEgenes <- (data05$padj<=0.05)
	genesDE<-data05[keepAllDEgenes,]
	print (c("Number of all DE genes:", dim(genesDE)[1]))
	write.table(genesDE, paste0("All.DE.genes.",conditionsName[2],"_vs_",conditionsName[(i)],".tab"), sep="\t", quote=F)
	#Obtainig genes with logFoldchage 1.5
	##keepAllDEgenesLFC2 <- (data05$padj<=0.05 & abs(data05$log2FoldChange) >= 1.5)
	##genesDELFC2<-data05[keepAllDEgenesLFC2,]
	#new way for filtering using dplyr
	DF05 <- as.data.frame(res05)
	DF05$Gene_id <- row.names(res05)
	DF05 <- DF05[,c(7,1:6)]
	genesDELFC1.5 <- dplyr::filter(DF05, padj <= 0.05 & abs(log2FoldChange) >= 0.585)
	print (c("Number of DE genes with LFC_1.5:", dim(genesDELFC1.5)[1]))
	write.table(genesDELFC1.5, paste0("All.DE.genes.LFC1.5.",conditionsName[2],"_vs_",conditionsName[(i)],".tab"), sep="\t", quote=F,row.names=F)
	#Generating MA plots
	pdf(paste0("MA_plot.",conditionsName[2],"_vs_",conditionsName[(i)],".pdf"))	
	plotMA(results(ddi, contrast=c("treatment_tissue",conditionsName[2],conditionsName[(i)]),alpha=0.05), main=paste0("DE Genes from:", conditionsName[1],"_vs_",conditionsName[(i)]))
	abline(h=c(-1.5,1.5),col="dodgerblue",lwd=2)	
	dev.off()
volcan <- EnhancedVolcano(res05,lab = "",x = "log2FoldChange",y = "pvalue", FCcutoff = 1.5,transcriptPointSize = 1.5,transcriptLabSize = 3.0,
        title = paste0("DE Genes from:", conditionsName[1],"_vs_",conditionsName[(i)]))
	ggsave(volcan, file=paste0("volcano.plot.",conditionsName[2],"_vs_",conditionsName[(i)],".png"))
}
##Start GNM WTM
print (c(conditionsName[4],"vs",conditionsName[6]))
#As we will use a p-adjust value cut-off <= 0.05 DESeq authors recomend an indepent filtering alpha=0.05 (https://www.pnas.org/content/107/21/9546.long)
res05 <- results(ddi, contrast=c("treatment_tissue",conditionsName[4],conditionsName[6]),alpha=0.05)
res05Ordered <- res05[order(res05$pvalue),]
data05 <- na.omit(as.data.frame(res05Ordered))
#writing all results
write.table(res05Ordered, paste0("DESEQ.",conditionsName[4],"_vs_",conditionsName[6],".tab"), sep="\t", quote=F)
print(summary(res05))
#Getting differential expressed genes below p-adjust value 0.05
keepAllDEgenes <- (data05$padj<=0.05)
genesDE<-data05[keepAllDEgenes,]
print (c("Number of all DE genes:", dim(genesDE)[1]))
write.table(genesDE, paste0("All.DE.genes.",conditionsName[4],"_vs_",conditionsName[6],".tab"), sep="\t", quote=F)
#Obtainig genes with logFoldchage 1.5
##keepAllDEgenesLFC2 <- (data05$padj<=0.05 & abs(data05$log2FoldChange) >= 1.5)
##genesDELFC2<-data05[keepAllDEgenesLFC2,]
#new way for filtering using dplyr
DF05 <- as.data.frame(res05)
DF05$Gene_id <- row.names(res05)
DF05 <- DF05[,c(7,1:6)]
genesDELFC1.5 <- dplyr::filter(DF05, padj <= 0.05 & abs(log2FoldChange) >= 0.585)
print (c("Number of DE genes with LFC_1.5:", dim(genesDELFC1.5)[1]))
write.table(genesDELFC1.5, paste0("All.DE.genes.LFC1.5.",conditionsName[4],"_vs_",conditionsName[6],".tab"), sep="\t", quote=F,row.names=F)
#Generating MA plots
pdf(paste0("MA_plot.",conditionsName[4],"_vs_",conditionsName[6],".pdf"))	
plotMA(results(ddi, contrast=c("treatment_tissue",conditionsName[4],conditionsName[6]),alpha=0.05), main=paste0("DE Genes from:", conditionsName[4],"_vs_",conditionsName[6]))
abline(h=c(-1.5,1.5),col="dodgerblue",lwd=2)	
dev.off()
volcan <- EnhancedVolcano(res05,lab = "",x = "log2FoldChange",y = "pvalue", FCcutoff = 1.5,transcriptPointSize = 1.5,transcriptLabSize = 3.0,
        title = paste0("DE Genes from:", conditionsName[1],"_vs_",conditionsName[(i)]))
	ggsave(volcan, file=paste0("volcano.plot.",conditionsName[4],"_vs_",conditionsName[6],".png"))
print ("Finishing Midgut comparisons")
#for PCA using VST form dds object generated above
vsd <- vst(ddi, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("treatment","tissue"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
PCA_plot <- ggplot(pcaData, aes(PC1, PC2, color=treatment,shape=tissue, label=pcaData$name)) + geom_point(size=3.5)+ geom_text(size=3)+xlab(paste0("PC1: ",percentVar[1],"% variance"))+ylab(paste0("PC2: ",percentVar[2],"% variance"))+coord_fixed()+scale_y_continuous(limits=c(-30,30), breaks=seq(-30,30,20))
ggsave(PCA_plot, file="PCA.pdf")
