"""
usage:
    python combine_sfs.py s1.csv s2.csv , ..., outputfile

Script to aggregate (SFS) the SFS.

Sums across the regions of the genome (by sample_id)
    Then summarize across replicates (samples: random_seed),
    using summary stats
"""
import sys

import numpy as np
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

# Summarize the replicates (samples from the trees)
# with quantiles


def quantile_q(quant):
    """
    Quantile factory function
    Args:
        quant: a number bewtween 0 and 1
    """
    return lambda x: np.quantile(x, quant)


def format_quantile_name(quant):
    """
    Helper function to construct column name
    Args:
        quant: a number bewtween 0 and 1
    """
    return 'Freq_q_' + str(quant).split('.')[1]


quantiles = [0.05, 0.25, 0.5, 0.75, 0.95]
quantiles = {format_quantile_name(q): quantile_q(q) for q in quantiles}
quantiles['min'] = np.min
quantiles['max'] = np.max

grp_vars.remove('random_seed')


(
    agg_data
    .groupby(grp_vars, as_index=False)
    .Freq
    .aggregate(quantiles)
    .to_csv(output, index=False)
)
