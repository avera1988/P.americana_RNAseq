#!/usr/bin/Rscript
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
load("GO.ancestry.RData")
args <- commandArgs(TRUE)
ANOT.total <- args[1]
OUT.total <- args[2]
UP <- args[3]
UP.out <- args[4]
DOWN <- args[5]
DOWN.out <- args[6]
annot.total <- read.delim(ANOT.total,sep="\t",header=T)
new_annot <- annot.total %>% 
  dplyr::left_join(GO, by="Gene_id")
write.csv(new_annot,OUT.total,row.names = F)
UP.total <- read.delim(UP,sep="\t",header=T)
new_up <- UP.total %>% 
  dplyr::left_join(GO, by="Gene_id")
write.csv(new_up,UP.out,row.names=F)
DOWN.total <- read.delim(DOWN,sep="\t",header=T)
new_down <- DOWN.total %>% 
  dplyr::left_join(GO, by="Gene_id")
write.csv(new_down,DOWN.out,row.names=F)
