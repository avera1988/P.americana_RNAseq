#!/usr/bin/Rscript
####################################
# This scripts produces VennDiagrams from All DE annotated Genes or Isoforms
# Author Aruro Vera 
# It needs the final annotation tables from Annoted.DE.v.6.sh
#
#################################################################


library(systemPipeR)
library(tidyr)
library(dplyr)
load("/home/avera/scripts/Annotation_RNASeq/Annotation.gene.isoform.cazy.AMP.Table.RData")
load("/home/avera/scripts/Analysis_pathways/DEseqTotal.isoforms.all.RData")
colnames(DEseqTotal)[1] <- "transcript_id"
row.names(DEseqTotal) <- DEseqTotal$transcript_id
print ("Starting Pipeline hindgut")
list.filenames <- list.files(pattern = "hindgut.tab.annot.final.tab$")
list.data<-list()
for (i in 1:length(list.filenames)){
  list.data[[i]]<-read.table(list.filenames[i],sep="\t",header=T,quote="")
}
names(list.data) <- list.filenames
list.data2 <- lapply(list.data,function(x){
  dplyr::select(x,contains("Isoform_id"))
  })
FNh=as.vector(unlist(list.data2[[1]]))
FWh=as.vector(unlist(list.data2[[2]]))
NWh=as.vector(unlist(list.data2[[3]]))
setlistH <- list(FN=FNh,FW=FWh,NW=NWh)
vennsetH <- overLapper(setlistH[1:3],type="vennsets")
pdf("Venn.Diagram.Hindgut.comparison.pdf",width = 9,height = 9)
vennPlot(vennsetH, mymain="Sharing DE genes in P. americana Hindgut")
dev.off()

print ("Starting Pipeline midgut")
list.filenames <- list.files(pattern = "midgut.tab.annot.final.tab$")
list.data<-list()
for (i in 1:length(list.filenames)){
  list.data[[i]]<-read.table(list.filenames[i],sep="\t",header=T,quote="")
}
names(list.data) <- list.filenames
list.data2 <- lapply(list.data,function(x){
  dplyr::select(x,contains("Isoform_id"))
})
FNm=as.vector(unlist(list.data2[[1]]))
FWm=as.vector(unlist(list.data2[[2]]))
NWm=as.vector(unlist(list.data2[[3]]))
setlistM <- list(FN=FNm,FW=FWm,NW=NWm)
vennsetM <- overLapper(setlistM[1:3],type="vennsets")
pdf("Venn.Diagram.Middgut.comparison.pdf",width = 9,height = 9)
vennPlot(vennsetM, mymain="Smaring DE genes in P. americana Midgut")
dev.off()

Unionmidgut <- as.data.frame(Unionmidgut)
colnames(Unionmidgut) <- "transcript_id"
Unionhindgut <- as.data.frame(Unionhindgut)
colnames(Unionhindgut) <- "transcript_id"

UnionmidgutDEG <- dplyr::inner_join(Unionmidgut,DEseqTotal,by="transcript_id") %>%
  dplyr::select(matches("transcript_id|FNM|FWM|NWM"))
UnionhindgutDEG <- dplyr::inner_join(Unionhindgut,DEseqTotal,by="transcript_id") %>%
  dplyr::select(matches("transcript_id|FNH|FWH|NWH"))

UnionmidgutDEGAnnot <- dplyr::inner_join(UnionmidgutDEG,Trinotate,by="transcript_id")
  

UnionhindgutDEGAnnot <- dplyr::inner_join(UnionhindgutDEG,Trinotate,by="transcript_id")


UnionMidgutHindgutBacteria <- list(UnionmidgutDEGAnnot,UnionhindgutDEGAnnot)
names(UnionMidgutHindgutBacteria) <- c("DE.Isoforms.Migut","DE.Isoforms.Hindgut")

openxlsx::write.xlsx(UnionMidgutHindgutBacteria,file="Bacterial.DE.isoforms.xlsx")
