"""
This scripts does the following:

    1. Extracts the data for a particular individual
    2. Collapses windows into continuous tracts: see -> https://github.com/santiago1234/mxb-genomes/blob/main/mxbgenomes/localancestry/tractslength.py
    3. save data.

usage:
    python scripts/process-individuals.py <mspfile> <outdir> <samplename>
    # this is an example, change the parameters accordingly
    python scripts/process-individuals.py data/3-pops/all-msp.csv data/3-pops/bed/NA19795.bed
"""

import pandas as pd
import os
import sys
sys.path.append('../../')
from mxbgenomes.localancestry import processxgmixout, tractslength
from mxbgenomes.utils import load_populations_info



# The msp file that contains all the samples
msp_file = sys.argv[1]
outdir = sys.argv[2]
individual = sys.argv[3]
msp = pd.read_csv(msp_file)


print('extracting: {}'.format(individual))
ind_d = processxgmixout.get_individual(msp, individual)
ind_d = processxgmixout.tidy_individual(ind_d)
# save uncollapsed
collapased_tracts = (
        ind_d
        .groupby(['chrn', 'Haplotype'])
        .apply(tractslength.collapse_windows_to_tracts)
        .reset_index(drop=True)
        .drop(['variable'], axis=1)
    )
outfile = os.path.join(outdir, individual + '.bed')
collapased_tracts.to_csv(outfile, sep='\t', index=False)
