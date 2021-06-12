#!/usr/bin/RScript
##################################################################
# This script merge all tables with Differential expression values
# after DESEq analysis
# Author Arturo Vear
# Sep 2019
#
#################################################################


library(dplyr)
library(tidyr)
library(tibble)
#Midgut DESEQ Table
DESEQFNM <- read.delim(file = "DESEQ.F_midgut_vs_N_midgut.tab", 
                       fill= T, h = T, sep = "\t")
DESEQFWM <- read.delim(file = "DESEQ.F_midgut_vs_WT_midgut.tab", 
                       fill= T, h = T, sep = "\t")
DESEQNWM <-read.delim(file = "DESEQ.N_midgut_vs_WT_midgut.tab", 
                      fill= T, h = T, sep = "\t")
DESEQFNM$Gene_id <- row.names(DESEQFNM)
DESEQFWM$Gene_id <- row.names(DESEQFWM)
DESEQNWM$Gene_id <- row.names(DESEQNWM)
DEseqM <- dplyr::full_join(DESEQFNM[,c(2,5,6,7)],DESEQFWM[,c(2,5,6,7)], by = "Gene_id",
                           suffix = c(".FNM", ".FWM")) %>%
  dplyr::full_join(DESEQNWM[,c(2,5,6,7)],by = "Gene_id") %>%
  dplyr::rename(padj.NWM = padj) %>% 
  dplyr::rename(log2FoldChange.NWM = log2FoldChange)%>%
  dplyr::rename(pvalue.NWM=pvalue) %>%
  dplyr::mutate_at(vars(log2FoldChange.FNM,log2FoldChange.FWM,log2FoldChange.NWM),
                   funs(round(., 2)))
DEseqM <- DEseqM[,c(4,1,2,3,5,6,7,8,9,10)]
#save(DEseqM,file="DESEQ.Midgut.total.RData")
rm(DESEQFNM,DESEQFWM,DESEQNWM)
#Hidgut DESEQ table
DESEQFNH <- read.delim(file = "DESEQ.F_hindgut_vs_N_hindgut.tab", 
                       fill= T, h = T, sep = "\t")
DESEQFWH <- read.delim(file = "DESEQ.F_hindgut_vs_WT_hindgut.tab", 
                       fill= T, h = T, sep = "\t")
DESEQNWH <-read.delim(file = "DESEQ.N_hindgut_vs_WT_hindgut.tab", 
                      fill= T, h = T, sep = "\t")
DESEQFHFM <- read.delim(file="DESEQ.F_hindgut_vs_F_midgut.tab", 
                        fill= T, h = T, sep = "\t")
DESEQNHNM <- read.delim(file="DESEQ.N_hindgut_vs_N_midgut.tab", 
                        fill= T, h = T, sep = "\t")
DESEQWHWM <- read.delim("DESEQ.WT_hindgut_vs_WT_midgut.tab", 
                        fill= T, h = T, sep = "\t")

DESEQFNH$Gene_id <- row.names(DESEQFNH)
DESEQFWH$Gene_id <- row.names(DESEQFWH)
DESEQNWH$Gene_id <- row.names(DESEQNWH)
DESEQFHFM$Gene_id <- row.names(DESEQFHFM)
DESEQNHNM$Gene_id <- row.names(DESEQNHNM)
DESEQWHWM$Gene_id <- row.names(DESEQWHWM)
DEseqH <- dplyr::full_join(DESEQFNH[,c(2,5,6,7)],DESEQFWH[,c(2,5,6,7)], by = "Gene_id",
                           suffix = c(".FNH", ".FWH")) %>%
  dplyr::full_join(DESEQNWH[,c(2,5,6,7)],by = "Gene_id") %>%
  dplyr::rename(padj.NWH = padj) %>% 
  dplyr::rename(log2FoldChange.NWH = log2FoldChange)%>%
  dplyr::rename(pvalue.NWH=pvalue) %>%
  dplyr::mutate_at(vars(log2FoldChange.FNH,log2FoldChange.FWH,log2FoldChange.NWH),
                   funs(round(., 2)))
DEseqH <- DEseqH[,c(4,1,2,3,5,6,7,8,9,10)]
rm(DESEQFNH,DESEQFWH,DESEQNWH)

DEseqBacteria <- dplyr::full_join(DEseqM,DEseqH,by="Gene_id")
rm(DEseqM, DEseqH)
DEseqTissues <- dplyr::full_join(DESEQFHFM[,c(2,5,6,7)],DESEQNHNM[,c(2,5,6,7)],by="Gene_id",
                                 suffix=c(".FHFM",".NHNM"))%>%
  dplyr::full_join(DESEQWHWM[,c(2,5,6,7)],by = "Gene_id") %>%
  dplyr::rename(padj.WHWM = padj) %>% 
  dplyr::rename(log2FoldChange.WHWM = log2FoldChange)%>%
  dplyr::rename(pvalue.WHWM=pvalue) %>%
  dplyr::mutate_at(vars(log2FoldChange.FHFM,log2FoldChange.NHNM,log2FoldChange.WHWM),
                   funs(round(., 2))) %>%
  dplyr::select(Gene_id,everything())
rm(DESEQFHFM,DESEQNHNM,DESEQWHWM)
DEseqTotal <- full_join(DEseqBacteria,DEseqTissues,by="Gene_id")
save(DEseqTotal,file="DEseqTotal.isoforms.all.RData")
