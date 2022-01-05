"""
Compute mL from countext counts.

usage:
	python compute-mL.py <context-counts> <mutation_rate_methylation_bins.txt> <outfile>
NOTE: CpG sites are excluded.
"""
import sys
import numpy as np
import pandas as pd


infile, ratesfile, outfile_ml = sys.argv[1:]

mL = pd.read_csv(ratesfile, sep='\t')
counts = pd.read_csv(infile, sep='\t')

def revcomplement(seq):
    pair = {
            'A': 'T', 
            'G': 'C',
            'T': 'A',
            'C': 'G'
    }
    return ''.join([pair[x] for x in seq[::-1]])

# We want to exclude the CpG sites
CpGs = [x + 'CG' for x in 'ACGT'] + ['CG' + x for x in 'ACGT']


mL_rev = mL.copy()
mL_rev['context'] = mL_rev.context.apply(revcomplement)
mL_rev['ref'] = mL_rev.ref.apply(revcomplement)
mL_rev['alt'] = mL_rev.alt.apply(revcomplement)

# drop the CpGs sites
mL = mL[~mL.context.isin(CpGs)]
mL_rev = mL_rev[~mL_rev.context.isin(CpGs)]


counts_rates = pd.concat([
    counts.merge(mL, on='context', how='inner'),
    counts.merge(mL_rev, on='context', how='inner')
    ])

counts_rates['mL'] = counts_rates.N * counts_rates.mu_snp

mr = counts_rates.mL.sum()

np.savetxt(outfile_ml, np.array([mr]))
