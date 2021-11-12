"""
Simulate random maiting population
Args:
    parent_tracts: dir containing the tracts data for the parental population.
    n_ind: the number of individuals to simulate in the offspring generation.
    seed: int, seed for random state.
    outdir: where to save simulted tracts data.
usage: python simulate.py <parent_tracts> <n_ind> <seed> <outdir>
"""

import sys
import os
from mating import utils, randmating

parent_tracts, n_ind, seed, outdir = sys.argv[1:]

if not os.path.exists(parent_tracts):
    sys.exit('dir: {} does not extis!'.format(parent_tracts))


tractspop = utils.load_tracts_pop(parent_tracts)

randmating.simulate_random_maiting_pop(
    tractspop, outdir, int(n_ind), int(seed))
