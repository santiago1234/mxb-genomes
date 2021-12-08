import pandas as pd

rates = pd.read_csv("data/mutation_rate_methylation_bins.txt", sep="\t")

def complement(seq):
    pair = {
            'A': 'T', 
            'G': 'C',
            'T': 'A',
            'C': 'G'
    }
    return ''.join([pair[x] for x in seq])


## PREPROCESSING

## ???????????????????????????
## TODO: check ???????????????
## get methylation_level = 0

rates = rates[rates.methylation_level == 0]

## ???????????????????????????
## TODO: check ???????????????
## total rate in triplet


rates = (
        rates.groupby(['context', 'ref'])
        ['mu_snp'].sum()
        .reset_index()
    )

rates_comp = rates.assign(context = lambda x: x['context'].apply(complement))


###

counts = pd.read_csv('data/counts/intergenic.csv')

counts = counts.rename(columns={'kmer': 'context'})


counts_rates = pd.concat([
    counts.merge(rates, on='context', how='inner'),
    counts.merge(rates_comp, on='context', how='inner')
    ])


counts_rates['total_mu'] = counts_rates['count'] * counts_rates['mu_snp']

print(counts_rates.total_mu.sum())
