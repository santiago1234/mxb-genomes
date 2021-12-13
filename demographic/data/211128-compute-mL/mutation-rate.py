""" Compute the mutation rate parameter

usage:
    python mutation-rate.py <mutation_rate_methylation_bins.txt> <counts>
"""
import pandas as pd
import sys

file_rate, file_count = sys.argv[1:]

rates = pd.read_csv(file_rate, sep="\t")
counts = pd.read_csv(file_count)

def revcomplement(seq):
    pair = {
            'A': 'T', 
            'G': 'C',
            'T': 'A',
            'C': 'G'
    }
    return ''.join([pair[x] for x in seq[::-1]])


rates = rates[rates.methylation_level == 0]

## total rate in triplet
rates = (
        rates.groupby(['context', 'ref'])
        ['mu_snp'].sum()
        .reset_index()
    )

rates_comp = rates.assign(context = lambda x: x['context'].apply(revcomplement))


# We want to exclude the CpG sites
CpGs = [x + 'CG' for x in 'ACGT'] + ['CG' + x for x in 'ACGT']

rates = rates[~rates.context.isin(CpGs)]
rates_comp = rates_comp[~rates_comp.context.isin(CpGs)]

### Compute mL

counts = counts.rename(columns={'kmer': 'context'})


counts_rates = pd.concat([
    counts.merge(rates, on='context', how='inner'),
    counts.merge(rates_comp, on='context', how='inner')
    ])


counts_rates['total_mu'] = counts_rates['count'] * counts_rates['mu_snp']

mL = counts_rates.total_mu.sum()

print(f'mL: {mL}')
