"""
Script to process simulation output.

This script will add the neutral mutations with a rate map using msprime.

NOTES:
    - Here, I take a random sample, pass a seed, of 49 individuals per population.
"""

import sys
import tskit
import os
import fwdpy11
import numpy as np
import pandas as pd
import itertools

sys.path.append('../../')
from simutils import utils, simulation

sim, out_file, path_to_samples, path_to_genetic_maps, graph, RANDOM_SEED = sys.argv[1:]

# Load data



RANDOM_SEED = int(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

# load the data used to set the simulation

# we can get the sample id from the simulation file
sample_id = os.path.basename(sim).split('-')[1]

simdat = utils.simuldata(path_to_samples, sample_id, path_to_genetic_maps)

ts = simulation.load_sim_as_ts(sim, graph)


# first we take a subsample from the three

N = 49
ts_s = simulation.subsample_individuals_pop(ts, N)


ts_by_mut_type = dict()
ts_by_mut_type['LOF'] = simulation.keep_selected_sites(ts_s, missense=False)
ts_by_mut_type['missense'] = simulation.keep_selected_sites(
    ts_s, missense=True)

# Add neutral mutations
ts_by_mut_type['noncoding'], ts_by_mut_type['synonymous'] = simulation.simulate_neutral_variation(
    ts_s, simdat)


def sfs(cat, pop):
    '''
    Get the folded SF for the given
    Mutation type: cat and Population: pop
    '''

    sf = simulation.get_single_pop_sfs(ts_by_mut_type[cat], pop)
    sf = sf.fold()

    d = {
        'Freq': sf.data,
        'Mask': sf.mask,
        'minor_allel_freq': range(sf.size)
    }

    sf = pd.DataFrame(d)

    sf = sf[~sf.Mask]
    sf['sample_id'] = sample_id
    sf['random_seed'] = RANDOM_SEED
    sf['Pop'] = pop
    sf['Mut_Type'] = cat

    return sf.loc[:, ['sample_id', 'random_seed', 'Pop', 'Mut_Type','minor_allel_freq', 'Freq']]


pops = ['YRI', 'IBS', 'CHB', 'MXB', 'MXL']
cats = ['missense', 'LOF', 'noncoding', 'synonymous']


data = pd.concat([sfs(x, y) for (x, y) in itertools.product(cats, pops)])
data.to_csv(out_file, index=False)
