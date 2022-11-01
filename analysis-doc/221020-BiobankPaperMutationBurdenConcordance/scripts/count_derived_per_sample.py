"""
Count the number of alternative alleles in
each individual of the passed vcf

usage:
    python count_derived_per_sample.py <VCF> <outfile.csv>
"""
import sys

import allel
import pandas as pd


def count_alternative(allele_counts):
    """
    Counts the number of alternative alleles
    Args:
        allele_counts: allel.AlleleCountsArray
    Returns: counts (int)
    """
    # 2nd column has the count for the alternative
    return allele_counts[:, 1].sum()


def main(callset, outfile):
    """
    Count the number of alternative alleles in
        each individual of the passed vcf
    Args:
        callset: VCF file
        outfile: output file with counts
    Returns:
        None, saves output file 
    """

    vcf = allel.read_vcf(callset)
    samples_to_index = {s: [i] for i, s in enumerate(vcf['samples'])}

    g_arr = allel.GenotypeArray(vcf['calldata/GT'])
    counts_per_sample = g_arr.count_alleles_subpops(samples_to_index,
                                                    max_allele=1)

    counts_per_sample = [(x, count_alternative(counts_per_sample[x]))
                         for x in counts_per_sample]

    counts_per_sample = pd.DataFrame(counts_per_sample,
                                     columns=['sample', 'alt_counts'])

    counts_per_sample.to_csv(outfile, index=False)


if __name__ == '__main__':
    par1, par2 = sys.argv[1:]
    main(par1, par2)
