"""
Compute the weighted mean for the b-scores.
The weights are the span of each interval in base pairs.

Usage:
    python weighted_mean.py <bfile>
"""
import pandas as pd
import sys

bfile = sys.argv[1]

# extract the chunk number
chunk = bfile.split('chunk_')[1].split('.')[0]

bfile = pd.read_csv(bfile, sep='\t', header=None)
bfile.columns = ['chr', 'start', 'end', 'name', 'bscore']

# the weights is the span of each interval in base pairs
bfile['weight'] = bfile['end'] - bfile['start']

# compute the weighted mean
weighted_mean = (bfile['bscore'] * bfile['weight']).sum() / bfile['weight'].sum()
print(f'chunk_{chunk} {weighted_mean}')
