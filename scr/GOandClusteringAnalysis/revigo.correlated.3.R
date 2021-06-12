#!/usr/bin/Rscript
###################################################################################################
#  Script for parsing revigo tables and correlate their Gene_ID with GO term collapsed
# Author Arturo Vera. Jan, 2019
# Usage: Rscript revigo.correlated.2.R REVIGO.csv Union.Goseq.enriched
#
###################################################################################################
library(tidyverse)
args <- commandArgs(TRUE)
REVIGO <- args[1]
UNION <- args[2]
revigo <- read_csv(REVIGO)
revigoMod <- dplyr::select(revigo, matches("TermID|Name|Eliminated")) %>%
  mutate(Eliminated=gsub(" ","",Eliminated))
enriched <- read_delim(UNION,delim = "\t")
enriched$gene_ids <- gsub(" ","",enriched$gene_ids)
enriched <- dplyr::rename(enriched, TermID=category)
correlated <- dplyr::inner_join(revigoMod, enriched,by="TermID") %>% 
  dplyr::select(matches("TermID|Name|Eliminated|gene_ids"))
write_tsv(correlated, 
            paste0(REVIGO,".correlated.tab"))
