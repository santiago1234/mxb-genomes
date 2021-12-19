"""
Remove CpG sites from VCF.
"""

from Bio import SeqIO
from cyvcf2 import VCF, Writer

# ref_genome = "/data/users/smedina/data-resources/genomes/human_g1k_v37.fasta"
ref_genome = 'chr22.fasta'

CpGsites = [
        'ACG', 'CGT',
        'CCG', 'CGG',
        'GCG', 'CGC',
        'TCG', 'CGA'
]

genome = SeqIO.to_dict(SeqIO.parse(ref_genome, "fasta"))


# Let's print the focal SNP and the REF SNP

vcf_file = '../results/data/210713-HardyW-filters/1TGP_and_50MXB-chr22-snps-vep-mask-HW-GRCh38.vcf.gz'

i = 0
for variant in VCF(vcf_file):
    i += 1
    if i < 20:
        break
    ##
    print(f'chr = {variant.CHROM}')
