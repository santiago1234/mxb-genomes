"""
Make a file with for each population a list of samples
"""
import sys
import pandas as pd

sys.path.append('../../../')

from mxbgenomes import utils

popinfo = utils.load_populations_info('../../../')

POPS = ['YRI', 'IBS', 'CHB', 'MXB', 'MXL']

for pop in POPS:
    samples = popinfo[popinfo['Subpopulation'] == pop]['Samplename'].tolist()
    with open(f'data/{pop}_samples.txt', 'w') as f:
        f.write('\n'.join(samples))
