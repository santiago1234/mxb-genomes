import pandas as pd
import numpy as np
import itertools
from tractsmodels import utils


# *****************************************************
# Get fits to data for bootstrap 0, that means the data.
# *****************************************************

pops = ['CLM', 'MXL', 'PEL', 'PUR']
mdls = ['ppx_xxp', 'ppx_xxp_pxx', 'ccx_xxp', 'ppx_ccx_xxp', 'ppc']
path_to_tracts_output = 'results/inference/'
bootstrap = 0  # 0 means it is the data and not a bootstrap replicate

params = itertools.product([path_to_tracts_output], mdls, pops, [0])
params = list(params)


fits_and_data = pd.concat([utils.ancestry_data_with_fits(*x) for x in params])
fits_and_data.to_csv('results/fits-data.csv', index=False)


# *****************************************************
# 4pops models
# *****************************************************

mdls = ['ccxx_xxpp', 'ppxx_ccxx_xxpp']
params = itertools.product(['results/inference-MXL-4pops/'], mdls, [0])

fits_and_data = pd.concat([utils.ancestry_data_with_fits_4pops(*x) for x in params])

fits_and_data.to_csv('results/fits-data-4pops-MXL.csv', index=False)


# *********************************************************
# Draw ancestry fractions over time for best fitting models
# *********************************************************

# Which is the best model for each population?
best_mdls = {
    'CLM': ('results/inference/CLM-ppx_ccx_xxp-boot0_mig', 'ppx_ccx_xxp'),
    'PEL': ('results/inference/PEL-ppx_ccx_xxp-boot0_mig', 'ppx_ccx_xxp'),
    'PUR': ('results/inference/PUR-ppp_pxp-boot0_mig', 'ppp_pxp'),
    'MXL': ('results/inference-MXL-4pops/ppxx_ccxx_xxpp-boot0_mig', 'ppxx_ccxx_xxpp')

}

default_labs = ['NAT', 'EUR', 'AFR']
pop_labs = {
    'CLM': default_labs,
    'PEL': default_labs,
    'PUR': ['EUR', 'NAT', 'AFR'],
    'MXL': ['EUR', 'NAT', 'AFR'] + ['EAS']
}


def anc_props(pop):
    migmat = np.loadtxt(best_mdls[pop][0])
    ap = utils.ancestry_prop_over_time(migmat, pop_labs[pop])
    ap['Population'] = pop
    ap['Model'] = best_mdls[pop][1]
    # get the pulses
    pulses = utils.migration_pulses(migmat, pop_labs[pop])
    pulses['Population'] = pop
    pulses['Model'] = best_mdls[pop][1]
    return ap, pulses


d = [anc_props(x) for x in best_mdls.keys()]

anps = pd.concat([x[0] for x in d])
pulses = pd.concat([x[1] for x in d])
anps.to_csv('results/anc-props-best-mdls.csv', index=False)
pulses.to_csv('results/pulses-best-mdls.csv', index=False)
