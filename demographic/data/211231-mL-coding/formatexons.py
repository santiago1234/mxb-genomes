"""
Parse data to use with VEP
"""
from Bio import SeqIO

NUCS = set('ACGT')
gfile = '/Users/santiagomedina/tmp/GRCh38-chr22.fasta'
genome = SeqIO.to_dict(SeqIO.parse(gfile, 'fasta'))


def get_seq(genome, chrom, start_pos, end_pos):
    """Get the reference sequence"""
    ref = genome[chrom][start_pos: end_pos]
    return str(ref.seq)


def variants(ref):
    variantes = NUCS - set(ref)
    return [f'{ref}/{x}' for x in variantes]


def process_exon(exon, genome):

    chrom, start, end, *_ = exon.strip().split('\t')

    chrom = chrom.replace('chr', '')
    start = int(start)
    end = int(end)
    exon_region = range(start, end + 1)

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


oufile = open('data/exon-snps.txt', 'w')


fexons = 'data/regions/exones.bed'
fexons = open(fexons, 'r')

for exon in fexons.readlines():
    vars_in_exon = process_exon(exon, genome)
    oufile.writelines(vars_in_exon)

oufile.close()
