"""
Script to generate the poplabels to estimate the genome wide genealogies.
I will use the following populations:
    - MXL
    - MXB
    - IBS
    - CHB
    - YRI
"""
import numpy as np
import pandas as pd
import sys
sys.path.append('../../')
from mxbgenomes.utils import load_populations_info

pops = load_populations_info('../../')

pops_to_use = ['MXL', 'MXB', 'IBS', 'YRI', 'CHB']
new_colnames = {
        'Subpopulation': 'populations',
        'Samplename': 'sample',
        'Superpopulation': 'group'
    }

pops = (
        pops[pops.Subpopulation.isin(pops_to_use)]
        .rename(columns=new_colnames)
        .assign(sex=lambda x: np.nan)
        .drop(columns=['Gender'])
    )


pops.to_csv("data/poplabels.txt", index=False, sep=' ')
