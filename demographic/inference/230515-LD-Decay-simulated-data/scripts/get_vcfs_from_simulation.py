"""
This script extracts VCFs with neutral sites from simulation for each population and saves the recombination map.
Usage:
    python get_vcfs_from_simulation.py <sim_id>
"""

import os
import sys

import pandas as pd
import tskit
import fwdpy11

sys.path.append('../../simulations/')
from simutils import utils, simulation


def main(sim_id):
    # Path to simulation results
    sim = f'../../simulations/220728-Simulation-DFE-Demography/results/simulations/sim-{sim_id}-pop.bin'

    graph = '../../simulations/220728-Simulation-DFE-Demography/ADMIXTURE-MXL.yml'

    ts = simulation.load_sim_as_ts(sim, graph)

    # Paths for simulation data
    path_to_samples = '../../data/220404-SimulationData/data/samples/'
    path_to_genetic_maps = '../../../resources/genetic-maps/'
    gmask = '../../../resources/genome-masks/20160622.allChr.mask.bed'

    simdata = utils.simuldata(path_to_samples, sim_id, path_to_genetic_maps)

    ts_nocd, _ = simulation.simulate_neutral_variation(ts, simdata)
    ts_nocd = simulation.filter_masked_sites(ts_nocd, simdata, gmask)

    # Subsample 100 individuals
    ts_nocd = simulation.subsample_individuals_pop(ts_nocd, 100)

    pops = ['YRI', 'IBS', 'CHB', 'MXB', 'MXL']

    for pop in pops:
        pop_ids = simulation.get_individuals_from_pops(ts_nocd, [pop])
        ts_pop = ts_nocd.simplify(pop_ids)

        # Write the VCF
        out_vcf = f'data/vcfs/sim{sim_id}_{pop}.vcf'
        vcf_string = ts_pop.as_vcf()
        with open(out_vcf, 'w') as f:
            f.write(vcf_string)

    # Save the recombination map
    pos = simdata.rmap.position

    recomb_map = pd.DataFrame({'pos': pos})
    recomb_map['cM'] = recomb_map.pos.apply(lambda x: simdata.rmap.get_cumulative_mass(x))
    recomb_map['chr'] = 1
    recomb_map[['pos', 'chr', 'cM']].to_csv(f'data/recomb_map/sim_{sim_id}.tsv', sep='\t', index=False)

    # save intervals
    intervals = simdata.noncoding_intervals
    # we set the chromosome to 1
    intervals['chro'] = 1
    intervals.to_csv(f'data/intervals/sim_{sim_id}.tsv', sep='\t', index=False, header=False)


if __name__ == "__main__":

    sim_id = sys.argv[1]
    main(sim_id)

