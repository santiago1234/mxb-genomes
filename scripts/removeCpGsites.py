"""
Remove CpG sites from VCF.

usage:
    python removeCpGsites.py <input.vcf> <refgenome.fasta> <output.vcf> <stats.txt> 

Important:
    - vcf file and reference genome must be in the same assembly
    - vcf chromosome names should match the names in the reference genome
"""
import sys
from Bio import SeqIO
from cyvcf2 import VCF, Writer


vcf_file, ref_genome, out_vcf, stats = sys.argv[1:]

genome = SeqIO.to_dict(SeqIO.parse(ref_genome, "fasta"))
vcf = VCF(vcf_file)

# create a new vcf Writer using the input vcf as a template.
w = Writer(out_vcf, vcf)

CpGsites = [
    'ACG', 'CGT',
    'CCG', 'CGG',
    'GCG', 'CGC',
    'TCG', 'CGA'
]


######################### FUNCTIONS #################################

def reverse_complement(seq):
    pairings = {
        'A': 'T',
        'G': 'C',
        'C': 'G',
        'T': 'A',
        'N': 'N'
    }
    return ''.join([pairings[x] for x in seq])[::-1]


def get_focal_snp(genome, variant):
    """
    returns *X* where X is the focal SNP and
    * are the sourrounding nucleotides.
    For examples, if the focal snp is A,
    the output could be CAG if C and G surround A.
    Args:
        genome: dict mapping chromosme names to sequences (SeqRecord)
    """
    focal_snp = genome[variant.CHROM][variant.start - 1: variant.start + 2].seq
    return focal_snp


def is_CpG(focal_snp):
    '''
    is the given focal mutation (or the reverse complement)
    a CpG site?
    '''
    if focal_snp in CpGsites:
        return True
    if reverse_complement(focal_snp) in CpGsites:
        return True
    return False


########################### DONE FUNCTIONS ###########################


n_CpGs = 0
for variant in vcf:
    ##
    focal_snp = get_focal_snp(genome, variant)

    # Sanity check
    if focal_snp[1:2] != variant.REF:
        raise ValueError(
            f'REF allel does not match REF-GENOME at position {variant.start}')
    if is_CpG(focal_snp):
        n_CpGs += 1
    if not is_CpG(focal_snp):
        w.write_record(variant)


CpGs_stats = f'CpG variants removed: {n_CpGs}\n'
stats_file = open(stats, 'w')
stats_file.write(CpGs_stats)
stats_file.close()


vcf.close()
w.close()
