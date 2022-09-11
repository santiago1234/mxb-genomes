"""
Code to compute Fst
"""

from itertools import combinations
import sys

from pybedtools import BedTool
import moments
import pandas as pd

sys.path.append('../')
from simutils import (simulation, utils)

POPS = ['MXL', 'YRI', 'IBS', 'CHB', 'MXB']
POP_PAIRS = list(combinations(POPS, 2))


def process_ts(ts_sim, simdat, gmask, samples_per_pop=50):
    """
    Preprocess ts_sim for computing
        the Fst (neutral variation)

    Preprocessin steps:
        1. Subsample indiviuals samples_per_pop 
        2. Add neutral variation genetic variation
        3. Remove masked regions

    Args:
        ts_sim: simulated tree sequence
        simdat: simulation data object
        gmask: Bedtool GRCh38 genome mask
        samples_per_pop: Samples per population

    Returns:
        preprocessed tree sequence
    """

    ts_preproc = simulation.subsample_individuals_pop(ts_sim, samples_per_pop)
    # we take the noncoding (intergenic and itronic)
    # and drop the synonymous with _
    ts_preproc, _ = simulation.simulate_neutral_variation(ts_preproc, simdat)

    ts_preproc = simulation.filter_masked_sites(ts_preproc, simdat, gmask)

    return ts_preproc


def pops_to_node_ids(ts_d):
    """
    Maps populations to node ids
    Returns:
        dict: pop_name -> list of node ids
    """
    pop_to_ids = {}

    for pop in POPS:
        node_ids = simulation.get_individuals_from_pops(ts_d, [pop])
        pop_to_ids[pop] = node_ids

    return pop_to_ids


def two_d_sfs(ts_d, fold=True):
    """
    Compute the 2d joint SFS.
    Returns:
        dict mapping pop pairs to moments.Spectrum
    """
    joint_sfs = {}
    pop_to_ids = pops_to_node_ids(ts_d)

    for pop1, pop2 in POP_PAIRS:
        pair_set = pop_to_ids[pop1], pop_to_ids[pop2]

        jsfs = ts_d.allele_frequency_spectrum(sample_sets=pair_set,
                                              polarised=True,
                                              span_normalise=False)

        jsfs = moments.Spectrum(jsfs, pop_ids=[pop1, pop2])
        if fold:
            jsfs = jsfs.fold()

        joint_sfs[(pop1, pop2)] = jsfs

    return joint_sfs


def fst(ts_d):
    """
    Compute the all pairwise Fst using the folded joint spectrum
    """

    jsfs = two_d_sfs(ts_d, fold=True)
    fst_s = []

    for pop_pair, spectrum in jsfs.items():
        fst_s.append([*pop_pair, spectrum.Fst()])

    return pd.DataFrame(fst_s, columns=['Pop1', 'Pop2', 'Fst'])


def main(cmd_arguments):
    """
    Script code
    """
    sim_id, path_to_project_root, out_fst = cmd_arguments

    sim_dat, sim_ts, d_files = utils.load_simulation_output(
        path_to_project_root, sim_id)

    gmask = BedTool(d_files['genome_masks'])

    ts_preproc = process_ts(sim_ts, sim_dat, gmask)

    fst_d = fst(ts_preproc)
    fst_d['Sim_id'] = sim_id
    fst_d.to_csv(out_fst, index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
