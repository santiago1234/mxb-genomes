library(GenomicRanges)
library(magrittr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(Biostrings)


# load genomi data --------------------------------------------------------

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
hs_genome <- BSgenome.Hsapiens.UCSC.hg38

# we want only the autosomes

autosomes <- paste0("chr", 1:22)
seqlevels(txdb) <- autosomes


# extract intronic regions ------------------------------------------------


intrones <- intronicParts(txdb)
intrones <- intrones[width(intrones) > 3]


export.bed(intrones, "data/regions/introns.bed")
intronic_seqs <- getSeq(hs_genome, intrones)
# add a fake name
names(intronic_seqs) <- paste0("intron_", 1:length(intronic_seqs))


writeXStringSet(intronic_seqs, "data/regions/introns.fasta", format = "fasta")

# get intergenic regions for H.  sapiens ----------------------------------
# See: Bioinformatic Data Skills, page: 319

chrom_grngs <- as(seqinfo(txdb), "GRanges")
collpased_tx <- reduce(transcripts(txdb))
strand(collpased_tx) <- "*"
intergenic <- setdiff(chrom_grngs, collpased_tx)
intergenic <- intergenic[width(intergenic) > 3]

export.bed(intergenic, "data/regions/intergenic.bed")

intergenic_seqs <- getSeq(hs_genome, intergenic)
names(intergenic_seqs) <- paste0("intergenic_", 1:length(intergenic_seqs))

writeXStringSet(intergenic_seqs, "data/regions/intergenic.fasta", format = "fasta")