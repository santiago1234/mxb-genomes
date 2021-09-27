"""
Generate table with populations for SFS

I will sample only 40 individual per population

"""
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


dims = {
        '1d': ['MXL'],
        '5d': ['IBS', 'YRI', 'CHB', 'MXB', 'MXL']
}


# sample individuals

pops = pops.groupby('Population').sample(n=40, random_state=42)

subset = lambda dim: pops[pops.Population.isin(dims[dim])]

# save pop-list

outdir = 'data/populations'
if not os.path.exists(outdir):
    os.makedirs(outdir)
# save
subset('1d').to_csv(os.path.join(outdir, 'pops-single.csv'), index=False)
subset('5d').to_csv(os.path.join(outdir, 'pops-five.csv'), index=False)
