"""
Compute pairwise fst.

usage:
    python scripts/fst_sim.py <ID> <path_root> <out_fst_table>

Args:
    ID: simulation id, a number
    path_root: path to project root mxb-genomes
    out_fst_table: output path for table with fst data
"""

from itertools import combinations
import sys

from pybedtools import BedTool
import moments
import pandas as pd

sys.path.append('../')
from simutils import (simulation, utils)

POPS = ['MXL', 'YRI', 'IBS', 'CHB', 'MXB']
POP_PAIRS = sorted(combinations(POPS, 2))


def process_ts(ts_sim, simdat, gmask, samples_per_pop=50):
    """
    Preprocess ts_sim for downstream analysis

    Preprocessin steps:
        1. Subsample indiviuals samples_per_pop 
        2. Get the trees for each functional category
        3. Remove masked regions

    Args:
        ts_sim: simulated tree sequence
        simdat: simulation data object
        gmask: Bedtool GRCh38 genome mask
        samples_per_pop: Samples per population

    Returns:
        (ts_nonsynonymous, ts_synonymous, ts_noncoding)
    """

    ts_preproc = simulation.subsample_individuals_pop(ts_sim, samples_per_pop)

    trees = dict()

    # this tree has the nonsynonymous variation
    trees['nonsynonymous'] = ts_preproc
    # we take the noncoding (intergenic and itronic)
    # and drop the synonymous with _
    ts_nocd, ts_syn = simulation.simulate_neutral_variation(ts_preproc, simdat)

    trees['synonymous'] = ts_syn
    trees['noncoding'] = ts_nocd

    def filter_mask(a_tree):
        '''
        remove masked regions
        '''
        return simulation.filter_masked_sites(a_tree, simdat, gmask)

    trees = {vcat: filter_mask(tree) for vcat, tree in trees.items()}
    return trees


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


def fst(ts_d, category=None):
    """
    Compute the pairwise Fst using the folded joint spectrum
    """

    jsfs = two_d_sfs(ts_d, fold=True)
    fst_s = []

    for pop_pair, spectrum in jsfs.items():
        fst_s.append([*pop_pair, spectrum.Fst()])

    fst_data = pd.DataFrame(fst_s, columns=['Pop1', 'Pop2', 'fst'])

    if category:
        fst_data['category'] = category

    return fst_data


def main(cmd_arguments):
    """
    Script code
    """
    sim_id, path_to_project_root, out_fst = cmd_arguments

    sim_dat, sim_ts, d_files = utils.load_simulation_output(
        path_to_project_root, sim_id)

    gmask = BedTool(d_files['genome_masks'])

    trees_by_vcat = process_ts(sim_ts, sim_dat, gmask)

    fst_data = [fst(trees_by_vcat[vcat], vcat) for vcat in trees_by_vcat]
    fst_data = pd.concat(fst_data)

    fst_data['sim_id'] = sim_id
    fst_data['pop_pair'] = fst_data.Pop1 + 'x' + fst_data.Pop2
    fst_data = fst_data.drop(columns=['Pop1', 'Pop2'])

    fst_data.to_csv(out_fst, index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
