library(dplyr)
annot <- read.delim("Trinotate.invertebrate.amp.kegg.cazy.tab", sep = "\t", h = T, fill = T) %>%
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

save(Annotations, file="Annotation.w.cazy.Table.RData")
