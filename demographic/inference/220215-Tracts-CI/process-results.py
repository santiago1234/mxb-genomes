import pandas as pd
import itertools
from tractsmodels import utils


# *****************************************************
# Get fits to data for bootstrap 0, that means the data.
# *****************************************************

pops = ['CLM', 'MXL', 'PEL', 'PUR']
mdls = ['ppx_xxp', 'ppx_xxp_pxx', 'ccx_xxp', 'ppx_ccx_xxp']
path_to_tracts_output = 'results/inference/'
bootstrap = 0  # 0 means it is the data and not a bootstrap replicate

params = itertools.product([path_to_tracts_output], mdls, pops, [0])

fits_and_data = pd.concat([utils.ancestry_data_with_fits(*x) for x in params])

fits_and_data.to_csv('results/fits-data.csv', index=False)
