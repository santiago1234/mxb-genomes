"""
Count derived alleles.

The sample name columns contain the number of alternative allele counts.
How do I polarize?
    - Keep SNPS where ANC == REF
    - Polarize when ALT == ANC
    - Remove the other varians:
        - lower case ancestral (e.g. a, c, g, t)
        - .
        - N
        - -
counts file has the format name: counts-var{VarClass}-chr{chrn}-type{Regions}.csv

Returns:
    - Table with Derived allele counts per sample (sum across variants)
    - I do this for all varints and the rare variants (MAF < 5%)
      The column VarFreq indicates this information.
"""
import pandas as pd
import sys

countsfile, outfile = sys.argv[1:]

# get info from filename
_, vtype, chrn, region = countsfile.split('-')
vtype = vtype.replace("var", "")
chrn = chrn.replace("chr", "")
region = region.replace("type", "").replace('.csv', "")

countsfile
d = pd.read_csv(countsfile)

non_sample_cols = ['ID', 'POS', 'REF', 'ALT', 'ANC', 'ALT_count']
# colnames of samples
sample_cols = [x for x in d.columns if x not in non_sample_cols]
# Polarize for ancestral allels
# Cases:
# case1:
# REF == ANC (uppercase)

case1 = d[d.REF == d.ANC]

# case2:
# polarize ALT == ANC

case2 = d[d.ALT == d.ANC]


def polarize(alt_count):
    """
    Polarize the number of allel counts to
    the ancestral allele
    """
    if alt_count == 2:
        return 0
    if alt_count == 0:
        return 2
    else:
        return 1


case2.loc[:, sample_cols] = case2.copy().loc[:, sample_cols].applymap(polarize)


polarized = pd.concat([case1, case2])
polarized['derived_count'] = polarized.loc[:, sample_cols].sum(axis=1)

total_alleles = len(sample_cols) * 2

polarized['DAF'] = polarized['derived_count'] / total_alleles


def count_derived_alleles(pol):
    dev_counts = pol.loc[:, sample_cols].sum()
    dev_counts = (
        dev_counts
        .reset_index()
        .rename(mapper={'index': 'Samplename', 0: 'derived_count'}, axis=1)
    )
    # add colinfo
    dev_counts['variant'] = vtype
    dev_counts['chr'] = chrn
    dev_counts['Region'] = region
    return dev_counts

# generate count for rare variants and all variants 

rare = polarized[polarized.DAF < 0.05]

counts_rare = count_derived_alleles(rare)
counts_rare['VarFreq'] = 'Rare (DAF <= 5%)'

counts_all = count_derived_alleles(polarized)
counts_all['VarFreq'] = 'Common (DAF <= 100%)'

counts = pd.concat([counts_rare, counts_all])
counts.to_csv(outfile, index=False)
