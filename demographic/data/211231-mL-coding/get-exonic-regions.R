library(GenomicRanges)
library(magrittr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(rtracklayer)
library(Biostrings)


# load genomi data --------------------------------------------------------

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
hs_genome <- BSgenome.Hsapiens.UCSC.hg38.masked

# we want only the autosomes

autosomes <- paste0("chr", 22:22)
seqlevels(txdb) <- autosomes


exones <- exonicParts(txdb) %>% 
  reduce()


export.bed(exones, "data/regions/exones.bed")
