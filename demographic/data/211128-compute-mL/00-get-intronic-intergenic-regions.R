library(GenomicRanges)
library(magrittr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)


# load genomi data --------------------------------------------------------

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
hs_genome <- BSgenome.Hsapiens.UCSC.hg38

# we want only the autosomes

autosomes <- paste0("chr", 1:22)
seqlevels(txdb) <- autosomes


# extract intronic regions ------------------------------------------------


intrones <- intronicParts(txdb)
intrones <- intrones[width(intrones) > 3]

intrones <- reduce(intrones)
names(intrones) <- paste0('intron_', 1:length(intrones))
export.bed(intrones, "data/regions/introns.bed")



# get intergenic regions for H.  sapiens ----------------------------------
# See: Bioinformatic Data Skills, page: 319

chrom_grngs <- as(seqinfo(txdb), "GRanges")
collpased_tx <- reduce(transcripts(txdb))
strand(collpased_tx) <- "*"

intergenic <- setdiff(chrom_grngs, collpased_tx)
intergenic <- intergenic[width(intergenic) > 3]
names(intergenic) <- paste0('intergenic_', 1:length(intergenic))

export.bed(intergenic, "data/regions/intergenic.bed")
