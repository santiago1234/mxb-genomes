"""
Remove CpG sites from VCF.
"""

from Bio import SeqIO
from cyvcf2 import VCF, Writer
import warnings


# NOTE: I am not using the correct reference genome.
# ref_genome = "/data/users/smedina/data-resources/genomes/human_g1k_v37.fasta"
ref_genome = 'GRCh38-chr22.fasta'

CpGsites = [
    'ACG', 'CGT',
    'CCG', 'CGG',
    'GCG', 'CGC',
    'TCG', 'CGA'
]

genome = SeqIO.to_dict(SeqIO.parse(ref_genome, "fasta"))


# Let's print the focal SNP and the REF SNP

vcf_file = '../results/data/210713-HardyW-filters/1TGP_and_50MXB-chr22-snps-vep-mask-HW-GRCh38.vcf.gz'

match_ref = 0
not_match_ref = 0
is_N_in_ref = 0
i = 0


for variant in VCF(vcf_file):
    i += 1
    if i > 100:
        break
    ##
    ref_genome = genome[variant.CHROM][variant.start]
    focal_snp = genome[variant.CHROM][variant.start - 1: variant.start + 2].seq
    ref_allel = variant.REF
    print(focal_snp, ref_allel)

    if ref_genome == 'N':
        is_N_in_ref += 1
    elif ref_genome != ref_allel:
        not_match_ref += 1
        warnings.warn(
            f'REF allel does not match REF-GENOME at position {variant.start}')
    else:
        match_ref += 1


stats = {
    'no_match_ref': not_match_ref,
    'match_ref': match_ref,
    'is_N_in_ref': is_N_in_ref
}
print(stats)
