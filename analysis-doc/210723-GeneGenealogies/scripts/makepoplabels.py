"""
Script to generate the poplabels to estimate the genome wide genealogies.
I will use the following populations:
    - MXL
    - MXB
    - IBS
    - CHB
    - YRI
    - PEL (samples representative of NAT ancestry)
"""
import numpy as np
import pandas as pd
import sys
sys.path.append('../../')
from mxbgenomes.utils import load_populations_info


pops = load_populations_info('../../')


nat_1tgp = np.loadtxt('../../resources/1TGP-samples-meta-data/native-american.txt', dtype=np.object_)

pops_to_use = ['MXB', 'MXL', 'IBS', 'YRI', 'CHB', 'PEL']
new_colnames = {
        'Subpopulation': 'population',
        'Samplename': 'sample',
        'Superpopulation': 'group'
    }

pops = (
        pops[pops.Subpopulation.isin(pops_to_use)]
        .rename(columns=new_colnames)
        .assign(sex=lambda x: np.nan)
        .drop(columns=['Gender'])
    )


# drop pel samples that are not NAT
pel = pops[pops.population == 'PEL']
pel_not_nat = pel['sample'][~pel['sample'].isin(nat_1tgp)].tolist()
pops = pops[~pops['sample'].isin(pel_not_nat)] 


# save samples
pops.to_csv("data/poplabels.txt", index=False, sep=' ')
