import pandas as pd
import numpy as np
import os
import sys
sys.path.append('../../../')
from mxbgenomes import utils


popinfoall = utils.load_populations_info('../../../')

main_pops = ['IBS', 'MXB', 'YRI', 'CHB', 'MXL']
pops = popinfoall[popinfoall.Subpopulation.isin(main_pops)]


pops = (
        pops
        .loc[:, ['Samplename', 'Subpopulation']]
        .rename(columns={'Subpopulation': 'Population'})
    )

pops = pops.groupby('Population').sample(n=50, random_state=42)

pops.to_csv('data/popsinfo.csv', index=False)

