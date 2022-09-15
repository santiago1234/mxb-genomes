"""
For each pair of populations writes-down a table containig
50 randomly choosen samples for each Population.

usage:
    python pairwise_poplist.py <output_dir_path>

Notes:
    output files will be name output_dir_path/pair_pop1Pop1-pop2Pop2.csv
"""

import sys
from itertools import combinations

import pandas as pd

sys.path.append('../../../../')
from mxbgenomes.utils import load_populations_info

OUTDIR = sys.argv[1]
POPS = ['MXL', 'YRI', 'IBS', 'CHB', 'MXB']

pop_pairs = combinations(POPS, 2)
popinfo = load_populations_info('../../../../')


def sample_per_pair(pop1: str, pop2: str) -> pd.DataFrame:
    """
    Fileter Subpopulations in (pop1, pop2)
    and take a sample of size 50 from each pop.
    Saves the output table.
    """

    samples_pair = (
        popinfo
        [popinfo.Subpopulation.isin((pop1, pop2))]
        .groupby(['Subpopulation'], as_index=False)
        .apply(lambda g: g.sample(50))
        .rename(columns={'Subpopulation': 'Population'})
        .loc[:, ['Samplename', 'Population']]
    )

    out_name = f'{OUTDIR}/pair_pop1{pop1}-pop2{pop2}.csv'

    samples_pair.to_csv(out_name, index=False)

    return samples_pair


_ = [sample_per_pair(*pair) for pair in pop_pairs]
