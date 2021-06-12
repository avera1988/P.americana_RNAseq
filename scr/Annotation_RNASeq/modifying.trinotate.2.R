######################################################################
# This script modifies the Trinoate annotation table adding CAZY 
# and Insect blast values. It returns RData objects with Trinotate
# modified table.
#
# Author: Arturo Vera
# Dec 2018
######################################################################

library(tidyr)
library(dplyr)

#load Trinotate data
trinotate <- read.delim("trinotate.annotation.all.xls",sep="\t",header=T,quote="")
#load Cazy
cazy <- read.delim("longest_orfs.pep.CAZY.0.2.annotate.unique",sep="\t",header=F,quote="")
colnames(cazy) <- c("prot_id","CAZyME","CAZY_annotation")
#load Kegg
kegg <- read.delim("triniti.id.with.Kegg.KO.tab",sep="\t",header=F,quote = "")
colnames(kegg) <- c("transcript_id","KO_number","KEGG_annot")
#load Blastx Invertebrate
blastxInver <- read.delim("trinity.fa.blastx.invertebrate.out.annot",sep="\t",header=T,quote="")
names(blastxInver)[c(1,2,5)] <- c("transcript_id","InvertebrateID","Invertebrate_Annotation")
blastxInver <- blastxInver[,c(1,2,5)]
#Removing duplicates
blastxInver <- dplyr::distinct(blastxInver)
#load AMP blastx
AMP <- read.delim("Trinity.amp.blastx.out.60.mod.annotated.mod",sep="\t",quote="")
names(AMP)[c(1,2)] <- c("transcript_id","AMP")
AMP <- dplyr::distinct(AMP)
#merging
trinotateMerged <- dplyr::full_join(blastxInver,trinotate,by="transcript_id") %>% 
  dplyr::full_join(kegg,by="transcript_id") %>% 
  dplyr::full_join(cazy,by="prot_id") %>% 
  dplyr::full_join(AMP,by="transcript_id")

#Geting same order that original trinotate
#use for to know the column names for (i in 1:length(colnames(trinotateMerged))) {a <- colnames(trinotateMerged)[[i]]; print (c(a,i))}
#for (i in 1:length(colnames(trinotateFinal))) {a <- colnames(trinotateFinal)[[i]]; print(c(i,a))}
trinotateFinal <- trinotateMerged[,c(4,1,5:18,2,3,19:23)]
#Obtaining Intersting Fields
trinotateUsefull <- trinotateFinal[,c(1:3,8,17:23)]
trinotateUsefull <- rename(trinotateUsefull,Gene_id =X.gene_id)
#keeping BX nad PFAM usefuls
trinotateUsefulMod <- trinotateUsefull %>% 
  separate(sprot_Top_BLASTX_hit,sep="\\^",into = c("Uniref_BX_ID",NA,NA,NA,NA,"Uniref_Annot","TAX")) %>% 
  separate(TAX,sep=";",into="Uniref_Kingdom") %>% 
  separate(Pfam, sep="\\^",into=c("PFAM_ID",NA,"PFAM_Annot"))
#Substituyin NA by .
toto <- as.matrix(trinotateUsefulMod)
y <- which(is.na(toto)==TRUE)
toto[y] <- "."
rm(y)
toto2 <- as.data.frame(toto)
rm(toto)
trinotateUseTbl <- as.tbl(toto2)
rm(toto2)
#Removing those not annotated genes
trinoNoAnnotated <- trinotateUseTbl %>% dplyr::filter(Uniref_BX_ID ==".") %>%
  dplyr::filter(PFAM_ID ==".") %>% dplyr::filter(InvertebrateID ==".") %>%
  dplyr::filter(KO_number==".") %>% dplyr::filter(CAZyME == ".") %>%
  dplyr::filter(AMP == ".")
#Annotated genes
trinotateAnnotateComplete <- dplyr::setdiff(trinotateUseTbl, trinoNoAnnotated)
save(trinotateUsefull,trinotateAnnotateComplete,file="Trinotate.2.RData")
write.table(trinotateAnnotateComplete,"Trinotate.invertebrate.amp.kegg.cazy.2.tab",sep="\t",quote=F,row.names = F)
#Collapsing tables
annot <- read.delim("Trinotate.invertebrate.amp.kegg.cazy.2.tab", sep = "\t", h = T, fill = T) %>%
  dplyr::select(-transcript_id) %>%
  dplyr::distinct() %>%
  dplyr::group_by(Gene_id) %>%
  dplyr::mutate(Uniref_BX_ID = paste0(Uniref_BX_ID, collapse = "", sep = ";" )) %>% 
  dplyr::mutate(PFAM_ID = paste0(PFAM_ID, collapse = "", sep = ";" )) %>% 
  dplyr::mutate(InvertebrateID = paste0(InvertebrateID, collapse = "", sep = ";" )) %>% 
  dplyr::mutate(KEGG_ID=paste0(KO_number, collapse = "",sep=";")) %>%
  dplyr::mutate(CAZyME=paste0(CAZyME, collapse = "",sep=";")) %>%
  dplyr::mutate(Description = paste0(Uniref_Annot,PFAM_Annot,Invertebrate_Annotation,KEGG_annot,CAZyME,collapse = "", sep = ";" ))  %>%
  dplyr::ungroup() %>% 
  dplyr::select(Gene_id, Uniref_BX_ID,InvertebrateID, PFAM_ID, KEGG_ID,CAZyME,Description) %>%
  dplyr::distinct()
cleanup_ids <- function(id_string){
  newstring <- gsub(";+", ";",
                    gsub("\\.+", ".",
                         gsub("^[;.]+|[;.]+$", "", id_string))
  )
  return(newstring)
}
Annotations <- annot %>% 
  mutate(Uniref_BX_ID = cleanup_ids(Uniref_BX_ID)) %>%
  mutate(InvertebrateID = cleanup_ids(InvertebrateID)) %>%
  mutate(PFAM_ID = cleanup_ids(PFAM_ID)) %>% 
  mutate(KEGG_ID = cleanup_ids(KEGG_ID)) %>%
  mutate(CAZyME = cleanup_ids(CAZyME)) %>%
  mutate(Description = cleanup_ids(Description))

save(Annotations, file="Annotation.w.cazy.Table.2.RData")

