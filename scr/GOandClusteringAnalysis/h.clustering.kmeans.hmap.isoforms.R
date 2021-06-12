################
#  This script gives you a heatmap with clustering by k-means
# Author Arturo Vera
# March 2020
#
#################
#libraries
library(magrittr)
library(dplyr)
library(gridExtra)
library(grid)
library(ggpubr)
library(stringr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(functional)
library(gplots)
library(dendextend)
library(RColorBrewer)
#######################

#Annotation, DESEQ and normalized counts tables
load("Annotation.gene.isoform.cazy.AMP.Table.RData")
load("DDS.Normalized.counts.isoforms.RData")
#load("DEseqTotal.RData")
load("DEseqTotal.isoforms.all.RData")
TrinDes <- Trinotate %>%
  unite("Description",sep=",",
        c(Uniref_Annot,PFAM_Annot,Invertebrate_Annotation,KO_number,KEGG_annot,CAZyME,CAZY_annotation,AMP)) %>%
  select(transcript_id,Description)

TrinDes <- data.frame(lapply(TrinDes,
                             function(x){
                               gsub("\\.","",gsub("\\.,","",x))
                             }))
DDS <- as.data.frame(table_counts_normalized)
rm(table_counts_normalized)
DDS$transcript_id <- row.names(DDS)
colnames(DEseqTotal)[1] <- "transcript_id"
row.names(DEseqTotal) <- DEseqTotal$transcript_id
DEseqTotal <- DEseqTotal %>%
  select(matches("transcript_id|Annot|FNM|FWM|NWM|FNH|FWH|NWH"))
system("cat *.tab|cut -f 2|sort|uniq  > Union.DEG.iso.both")
tb <- read.delim(list.files(path = ".", pattern = "*both"),sep="\t",header = T)
tb <- tb %>%
  dplyr::rename("transcript_id"=Isoform_id) %>%
  mutate_if(is.factor,as.character)
tb <- as.data.frame(tb)
tbDDS <- dplyr::inner_join(tb,DDS,by="transcript_id")
row.names(tbDDS) <- tbDDS$transcript_id
tbDDSMat <- tbDDS[,2:19]
#PHTotal <- phMap(tbDDSMat)
MatrixLog2Mas1 <- log2(tbDDSMat+1)
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
data_z_score <- t(apply(MatrixLog2Mas1,1,cal_z_score))

#Clustering by k-means
set.seed(100)
Clustering <- pheatmap(data_z_score,kmeans_k = 5)
clusterDF <- as.data.frame(factor(Clustering$kmeans$cluster))
colnames(clusterDF) <- "Cluster"
OrderByCluster <- data_z_score[order(clusterDF$Cluster),]

ph_annot <- as.data.frame(colnames(OrderByCluster))
names(ph_annot) <- "names"
ph_annot$Tissue <- ifelse(grepl("H",ph_annot$names), "Hindgut","Midgut")
row.names(ph_annot) <- ph_annot$names
ph_annot$Bacterial_treatment <- ifelse(grepl("F",ph_annot$names),"Germfree",
                                       ifelse(grepl("N",ph_annot$names),"Gnotobiotic","Wildtype"))
ph_annot <- dplyr::select(ph_annot,Bacterial_treatment, Tissue)
Colors <- brewer.pal(5,"Set2")

annot_colors=list(Bacterial_treatment=c(Germfree = "#CCCCCC",
                                        Gnotobiotic = "#999999",
                                        Wildtype = "#333333"),
                  Tissue = c(Hindgut="#9933FF",Midgut="#3399FF"),
                  Cluster=c("1"="#66C2A5","2"="#FC8D62","3"="#8DA0CB","4"="#E78AC3","5"="#A6D854"))

pmap <- pheatmap(OrderByCluster,annotation_row = clusterDF,treeheight_col = 0,treeheight_row  =0,
                 show_rownames = FALSE,cluster_rows = FALSE,color=redgreen(n=190),fontsize_row=10,
                 cellwidth=18,annotation_col = ph_annot,annotation_colors = annot_colors)

ggsave(pmap, file="heatmap.DEG.total.isoforms.pdf",width = 12,height = 12,dpi=300)

#Getting clustering and annotation

clusterDF$transcript_id <- rownames(clusterDF)
clusterDFSorted <- dplyr::arrange(clusterDF, clusterDF$Cluster)
clusterDFSorted <- clusterDFSorted[,c(2,1)]

Nclusters=(Clustering$kmeans$iter)+1
#Mergin Clusters with Annotations

#This for generates n Objects as kmeans in the clusterign
for (i in 1:Nclusters){
  assign(paste("Gene_Cluster", i, sep = ""), subset(clusterDFSorted, clusterDFSorted$Cluster == i))
}
#To generate factors names of each cluster
Numbers <- 1:Nclusters
n <- length(Numbers)
geneC <- paste("Gene_Cluster", 1:n, sep="")
#This function returns each cluster into the different cluster
GetDatas <- function(x){
  dat <- get(x)
  #Rnames <- rownames(dat)
  #dat$Gene_id <- Rnames
  dat
}
geneDatas <- lapply(geneC, GetDatas)

#Getting annotations
GetAnnotation <- function(x){
  datos<- as.data.frame(x)
  mezcla <- dplyr::inner_join(datos,TrinDes,by="transcript_id")
  mezcla
}
clusterAnnotations <- lapply(geneDatas, GetAnnotation)
x <- ""
for(i in 1:n){ 
  x[i] <- paste0("Union.",i)
}

names(clusterAnnotations) <- x
list.Cluster.Annot.DEG <- lapply(clusterAnnotations,function(x){
  dplyr::inner_join(x,DEseqTotal,by="transcript_id")
})
Columns <- colnames(list.Cluster.Annot.DEG$Union.1)
#Adding LFC values using DEseq
list.Cluster.Annot.DEG <- lapply(clusterAnnotations,function(x){
  dplyr::inner_join(x,DEseqTotal,by="transcript_id") %>%
    distinct(transcript_id,.keep_all = T)
})
#Generating a vector with names
Columns <- colnames(list.Cluster.Annot.DEG$Union.1)

#Saving tables with all info for each cluster
for (i in 1:length(list.Cluster.Annot.DEG)){
  write.table(list.Cluster.Annot.DEG[i],
              file=paste0(names(list.Cluster.Annot.DEG)[i],
                          ".Cluster.Annotation.DEG.tab"),
              sep="\t",
              row.names = F,
              quote=F,
              col.names = Columns,
              fileEncoding="UTF-8")
}

