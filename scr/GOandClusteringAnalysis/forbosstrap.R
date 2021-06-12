###########################################
#
#  This script takes the z-score table and 
# perform a HC analysis and boostrap 
# 
# Dependencies: data_z_score from h.cluter.kmeans.R script
# Author: Arturo Vera
# May 2021
#
##########################################



set.seed(100)
d <- dist(data_z_score, method="euclidean") 

pfit <- hclust(d, method="complete")   
plot(pfit)

km = kmeans(data_z_score, 5, iter.max = 100)

library(fpc)

cboot.hclust <- clusterboot(data_z_score,clustermethod =kmeansCBI ,
                            k=5,seed = 100)
cboot.hclust$bootmean

clusterDF <- as.data.frame(factor(cboot.hclust$partition))
colnames(clusterDF) <- "Cluster"
OrderByCluster <- data_z_score[order(clusterDF$Cluster),]

ph_annot <- as.data.frame(colnames(OrderByCluster))
names(ph_annot) <- "names"
ph_annot$Tissue <- ifelse(grepl("H",ph_annot$names), "Hindgut",
                          "Midgut")
row.names(ph_annot) <- ph_annot$names
ph_annot$Bacterial_treatment <- ifelse(grepl("F",ph_annot$names),"Germfree",
                                       ifelse(grepl("N",ph_annot$names),
                                              "Gnotobiotic",
                                              "Wildtype"))
ph_annot <- dplyr::select(ph_annot,Bacterial_treatment, Tissue)
Colors <- brewer.pal(5,"Set2")

annot_colors=list(Bacterial_treatment=c(Germfree = "#CCCCCC",
                                        Gnotobiotic = "#999999",
                                        Wildtype = "#333333"),
                  Tissue = c(Hindgut="#9933FF",Midgut="#3399FF"),
                  Cluster=c("1"="#66C2A5",
                            "2"="#FC8D62",
                            "3"="#8DA0CB",
                            "4"="#E78AC3",
                            "5"="#A6D854"))

pmap <- pheatmap(OrderByCluster,
                 annotation_row = clusterDF,
                 treeheight_col = 0,
                 treeheight_row  =0,
                 show_rownames = FALSE,
                 cluster_rows = FALSE,
                 color=redgreen(n=190),
                 fontsize_row=10,
                 cellwidth=18,
                 annotation_col = ph_annot,
                 annotation_colors = annot_colors)
ggsave(pmap, 
       file="heatmap.DEG.total.isoforms.pdf",
       width = 8,
       height = 8,
       dpi=300)

#Getting clustering and annotation

clusterDF$transcript_id <- rownames(clusterDF)
clusterDFSorted <- dplyr::arrange(clusterDF, clusterDF$Cluster)
clusterDFSorted <- clusterDFSorted[,c(2,1)]

Nclusters=as.vector(levels(clusterDFSorted$Cluster))
#Mergin Clusters with Annotations

#This for generates n Objects as kmeans in the clusterign
for (i in 1:length(Nclusters)){
  assign(paste("Gene_Cluster", i, sep = ""), 
         subset(clusterDFSorted, 
                clusterDFSorted$Cluster == i))
}
#To generate factors names of each cluster
Numbers <- 1:length(Nclusters)
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
              file=paste0(names(list.Cluster.Annot.DEG)[i],".Cluster.Annotation.DEG.tab"),
              sep="\t",row.names = F,quote=F,col.names = Columns)
}
