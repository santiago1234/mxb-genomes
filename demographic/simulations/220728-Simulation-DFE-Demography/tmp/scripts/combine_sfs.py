"""
usage:
    python combine_sfs.py s1.csv s2.csv , ..., outputfile


Script to aggregate (SFS) the SFS.


Sums across the sample id, this to 
get the SFS that represents the simulated region

"""
import sys

import pandas as pd

*input_sfs, output = sys.argv[1:]


data = pd.concat([pd.read_csv(x) for x in input_sfs])

grp_vars = ['random_seed', 'Pop', 'Mut_Type', 'minor_allel_freq']


agg_data = (
    data
    .groupby(grp_vars)
    .Freq
    .sum()
    .reset_index()
)

agg_data.to_csv(output, index=False)
