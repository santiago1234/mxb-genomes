## See EnsDb doc: https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#4_Retrieving_sequences_for_genetranscriptexon_models
library(EnsDb.Hsapiens.v86)
library(magrittr)

edb <- EnsDb.Hsapiens.v86
edb <- 
  edb %>% 
  filter(SeqNameFilter(1:22))
exonsByGene <- exonsBy(edb, filter = GeneBiotypeFilter("protein_coding"), by = "gene")

## Create a GRanges of non-overlapping exon parts.
## Note that I am only retrieving 
DJE <- disjointExons(edb, filter = AnnotationFilterList(
  SeqNameFilter(c(1:22)),
  GeneBiotypeFilter("protein_coding")))

rtracklayer::export.bed(DJE, 'data/exons.bed')
