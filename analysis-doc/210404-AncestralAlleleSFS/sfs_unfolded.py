"""
Compute the unfolded SFS using the aa info.
Only valid for biallelic loci
"""

import allel
import pandas as pd
import numpy as np

aa_file = "ancestral.csv"
vcf_file = "variant_effect_output.vcf.gz"

vcf = allel.read_vcf(vcf_file, alt_number=1)
aa = pd.read_csv(aa_file)
aa = aa.drop_duplicates(subset=['ID'])

# Make a frame
df = {
    'ID': vcf['variants/ID'],
    'REF': vcf['variants/REF'],
    'ALT': vcf['variants/ALT']
}
df = pd.DataFrame(df)

# Add the ancestral allel

df = pd.merge(df, aa, how='left', on='ID')


df = df.set_index('ID')

# Sanity check, that the variants are aligned with the vcf

if not np.all(vcf['variants/ID'] == df.index):
    raise ValueError('sanity check not passed, data is no aligned')


# Possible scenarios
#     a) REF allele is the AA
#     b) neither REF or ALT are the AA
#     c) We do not know the Ancestral allel, this is indicated with NA or .
#     d) ALT is the AA
# For the cases a-c the REF will be considered the AA.
# Only in the case d) the REF will be switched with the AA
# to compute the unfolded SFS


# Conditions
is_alt_aa = (df.ALT == df.AA).to_numpy()

# Extract the genotype data

ga = vcf['calldata/GT']


ga_ref_is_aa = ga[~is_alt_aa, :, :]
ga_alt_is_aa = ga[is_alt_aa, :, :]

ga_alt_is_aa = allel.GenotypeArray(ga_alt_is_aa)
ga_ref_is_aa = allel.GenotypeArray(ga_ref_is_aa)

ac_alt_is_aa = ga_alt_is_aa.count_alleles(max_allele=1)
ac_ref_is_aa = ga_ref_is_aa.count_alleles(max_allele=1)

# Compute the SFS


def clear_variants(acounts, n):
    """
    Remove variants where:
        1) The sum of allele counts is not 2 * n. This
        happens because for some variants we could have
        missing genotypes.
        2) The variant is not a variant in the set of
        samples. In this case the counts for the ancestral
        allel is 0, or the counts for the derived alle is 0
    Args:
        acounts: allel counts array
        n: sample size
    """

    diploid_sample_size = 2 * n
    acounts = acounts[acounts.sum(axis=1) == diploid_sample_size]
    one_count_is_zero = np.logical_or(acounts[:, 0] == 0, acounts[:, 1] == 0)
    return acounts[np.logical_not(one_count_is_zero), ]


# clear variants
n_samples = vcf['samples'].size
ac_alt_is_aa = clear_variants(ac_alt_is_aa, n_samples)
ac_ref_is_aa = clear_variants(ac_ref_is_aa, n_samples)

# Compute SFS
sfs_ref_is_aa = allel.sfs(ac_ref_is_aa[:, 1])

# We pass the 
sfs_alt_is_aa = allel.sfs(ac_alt_is_aa[:, 0])

# Now we sum
sfs = sfs_ref_is_aa + sfs_alt_is_aa

