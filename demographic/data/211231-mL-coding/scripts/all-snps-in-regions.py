"""
Generate a bed file with all posible biallelic variants for each position in the
given bed file.
The output file is meant to be annotated with VEP.

usage:
    python scripts/all-snps-in-regions.py <regions.bed> <refgenome.fa> <all-variants-in-regions.txt>

Args:
    - regions.bed: Input bed file containing regions (for example exons).
    - refgenome.fa: The reference genome, same build as regions.bed. Also, chromosme names
        should not include the 'chr' prefix.
    - all-variants-in-regions.tx: file path for output data. The output has the VEP
        input format: https://www.ensembl.org/info/docs/tools/vep/vep_formats.html#default
"""
import sys
from Bio import SeqIO


fexons, gfile, outvariants = sys.argv[1:]


genome = SeqIO.to_dict(SeqIO.parse(gfile, 'fasta'))
NUCS = set('ACGT')


def get_seq(genome, chrom, start_pos, end_pos):
    """Get the reference sequence"""
    ref = genome[chrom][start_pos: end_pos]
    return str(ref.seq)


def variants(ref):
    """List of possible variant alleles"""
    variantes = NUCS - set(ref)
    return [f'{ref}/{x}' for x in variantes]


def process_exon(exon, genome):

    chrom, start, end, *_ = exon.strip().split('\t')

    chrom = chrom.replace('chr', '')
    start = int(start)
    end = int(end)
    exon_region = range(start - 1, end + 1)

    data = list()

    for pos in exon_region:
        ref_seq = get_seq(genome, chrom, pos, pos + 1)
        # get the context
        context = get_seq(genome, chrom, pos - 2, pos + 3)
        # make and id
        v_id = chrom + ':' + str(pos + 1) + ':' + '_' + context + '_:'
        d = [
            f'{chrom}\t{pos+1}\t{pos+1}\t{x}\t.\t{v_id}{x}\n' for x in variants(ref_seq)]
        data.extend(d)

    return data


oufile = open(outvariants, 'w')


fexons = open(fexons, 'r')

for exon in fexons.readlines():
    vars_in_exon = process_exon(exon, genome)
    oufile.writelines(vars_in_exon)

oufile.close()
