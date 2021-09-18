"""
Script to parse vcf 
usage:
    python processdata.py <VCF.vcf> <ANC.fa> <output.csv>
The VCF file should contain data for and and only one chromosome
The ANC genome file should only contain one and only one chromosome

Output:
    - A csv table with the following columns:
        - ID: variants id
        - POS: genomic position
        - REF: reference allele
        - ALT: alternative allele
        - ANC: ancestral allele
        - ALT_count: number of alternative allele counts per variant
        - SAMPLES (multiple columns): each entry records the number of alt alleles
        for each sample at the given variant (range 0-2).
"""
import pandas as pd
import allel
from Bio import SeqIO
import sys


vcf, anc_genome, outfile = sys.argv[1:]


# Load ancestral genome
anc_genome = SeqIO.read(anc_genome, 'fasta')
anc_seq = anc_genome.seq

vcf = allel.read_vcf(vcf, alt_number=1)
positions = vcf['variants/POS']
ga = allel.GenotypeArray(vcf['calldata/GT'])
samples = vcf['samples']
# map samples to position in calldata
samples_indices = {x:[i] for (i, x) in enumerate(samples)}


def count_alt_allels_per_sample_and_variant():
    # Returns a data frame with alternative allel counts
    # per sample and variants
    # get allele counts
    alt_counts = ga.count_alleles_subpops(subpops=samples_indices, max_allele=1)
    get_alt_count = lambda sample: alt_counts[sample][:, 1]
    counts = pd.DataFrame({s: get_alt_count(s) for s in samples})
    counts['ID'] = vcf['variants/ID']
    return counts


def get_variant_data():
    tmp = str(anc_seq)
    anc_allels = []
    for p in positions:
        anc_allels.append(tmp[p - 1])

    d = {
        'ID': vcf['variants/ID'],
        'POS': positions,
        'REF': vcf['variants/REF'],
        'ALT': vcf['variants/ALT'],
        'ANC': anc_allels,
        'ALT_count': ga.count_alleles(max_allele=1)[:, 1]
    }
    return pd.DataFrame(d)


d = pd.merge(
        get_variant_data(),
        count_alt_allels_per_sample_and_variant(),
        how='inner', on='ID')

d.to_csv(outfile, index=False)

