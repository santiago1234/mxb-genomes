## See EnsDb doc: https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#4_Retrieving_sequences_for_genetranscriptexon_models
library(EnsDb.Hsapiens.v86)
library(magrittr)

edb <- EnsDb.Hsapiens.v86
edb <- 
  edb %>% 
  filter(SeqNameFilter(1:22)) %>% 
  filter(GeneBiotypeFilter("protein_coding"))

cds <- cdsBy(edb)

cds %>% 
  unlist() %>% 
  reduce() %>% 
  rtracklayer::export.bed('data/coding-regions.bed')
