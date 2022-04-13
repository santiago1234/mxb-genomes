"""
Compute mL from vep output.
"""
import re
import pandas as pd
from Bio.Seq import Seq

vep_file = 'vep-chr22-UNIQ.txt.gz'
mus_file = '../211128-compute-mL/data/mutation_rate_methylation_bins.txt'
chrn = 22
output = 'mL-chr22.csv'

vep = pd.read_csv(vep_file, sep='\t')

mus = pd.read_csv(mus_file, sep='\t')

# NOTE: Filter data by methylation_level
mus = mus[mus.methylation_level == 0]

cpgs = ["CGA", "CGT", "CGG", "ACG", "TCG", "GCG", "CCG", "CGC"]
exclude_cpgs = True

# Consequences
# We only want to consider the following categories
Q = {
    'Q': ['missense', 'synonymous', 'LOF', 'LOF', 'LOF'],
    'Consequence': ['missense_variant', 'synonymous_variant', 'stop_lost', 'stop_gained', 'start_lost']
}

Q = pd.DataFrame(Q)

vep = pd.merge(vep, Q, on='Consequence', how='inner')

vep = vep.drop(columns=['Gene', 'Feature_type'])  #We dont need this info


# This function extracts the context sequence from the id
get_context = lambda x: re.search(r'_[ACGT]{5}_', x).group(0)[2:-2]

vep['context'] = vep.Uploaded_variation.map(get_context)


# Counts
counts = vep.groupby(['context', 'Allele', 'Q']).size().reset_index().rename(columns={0:'n'})

# For the mus table we want also the reverse complement

def rev_com(aseq):
    s = Seq(aseq)
    return str(s.reverse_complement())

mus_rev = mus.copy()

mus_rev['context'] = mus_rev.context.map(rev_com)
mus_rev['ref'] = mus_rev.ref.map(rev_com)
mus_rev['alt'] = mus_rev.alt.map(rev_com)

all_mus = pd.concat([mus, mus_rev])

if exclude_cpgs:
    all_mus =  all_mus[~all_mus.context.isin(cpgs)]


## Compute mL
counts = counts.rename(columns={'Allele': 'alt'})

mLs = pd.merge(counts, all_mus, on=['context', 'alt'], how='inner')
mLs['mL'] = mLs.mu_snp * mLs.n

mLs = mLs.groupby(['Q']).mL.sum().reset_index()
mLs.assign(chrn=chrn).to_csv(output, index=False)
