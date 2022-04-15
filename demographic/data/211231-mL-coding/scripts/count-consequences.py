"""
Script to count how many variants are synonymous, missense, and lof.
"""
import pandas as pd
import re
import sys

vep_file, chrn, output = sys.argv[1:]
vep = pd.read_csv(vep_file, sep='\t')

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


# exclude the cpgs
vep = vep[~vep.context.isin(cpgs)]

# now count
counts = vep.groupby(['Q']).size().reset_index().rename(columns={0:'n'})
counts.assign(chrn=chrn).to_csv(output, index=False)
