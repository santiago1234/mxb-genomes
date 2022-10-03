"""
Compute Fst scores in simulation output
usage:
    python scripts/fst_msprimesimul.py <TREES> <n_cores> <outfile>
"""
import os
import sys
import multiprocessing as mp

import pandas as pd
import tskit
import moments
import pybedtools

sys.path.append('scripts/')
sys.path.append('../220910-Fst-Comparison/scripts/')
sys.path.append('../')

from simutils import utils, simulation
from neutralsimulator import POPs
from fst_sim import POP_PAIRS


def map_pops_to_ids(tree):
    """
    Maps pop codes to individual ids
    It assumes the first K=50 come from POPs[0]
        second Kth POPs[1] 
        ...
    Returns:
        dict mapping pop codes to ids
    """
    n_pops = len(POPs)
    n_by_pop = tree.num_individuals // n_pops
    # since we are working with diploid data, we multiply
    # by a factor of 2
    n_by_pop *= 2
    inds_by_pop = {
        POPs[i]: range(i * n_by_pop, (i + 1) * n_by_pop)
        for i in range(n_pops)
    }
    inds_by_pop = {k: list(v) for k, v in inds_by_pop.items()}
    return inds_by_pop


def jsfs(tree):
    """
    Compute the joint SFS
    Args:
        tree: tskit.Tree
        popair: pop1, pop2
    Returns:
        joint sfs
    """
    spectrums = {}
    for popair in POP_PAIRS:
        pop1, pop2 = popair
        ids_by_pop = map_pops_to_ids(tree)
        pair_set = ids_by_pop[pop1], ids_by_pop[pop2]

        spectrum = tree.allele_frequency_spectrum(sample_sets=pair_set,
                                                  polarised=True,
                                                  span_normalise=False)

        spectrum = moments.Spectrum(spectrum, pop_ids=popair).fold()
        spectrums[popair] = spectrum

    return spectrums


def fst(tree, category=None):
    """
    Compute the pairwise Fst using the folded joint spectrum
    """

    spectrums = jsfs(tree)
    fst_s = []

    for pop_pair, spectrum in spectrums.items():
        fst_s.append([*pop_pair, spectrum.Fst()])

    fst_data = pd.DataFrame(fst_s, columns=['Pop1', 'Pop2', 'fst'])

    if category:
        fst_data['category'] = category

    return fst_data


def decode_metadata_from_filename(tree_file):
    """
    Retrieve metadata from file name
    """
    filename = os.path.basename(tree_file)
    simid = filename.split('sim_')[1].split('-cat')[0]
    vcat = filename.split('cat_')[1].split('.ts')[0]
    return simid, vcat


def get_fsts(tree_file, root_path='../../../'):
    """
    Args:
        root_path: relative path to project root
    """
    tree = tskit.load(tree_file)
    simid, vcat = decode_metadata_from_filename(tree_file)

    sim_data, _, d_files = utils.load_simulation_output(root_path, simid)
    gmask = pybedtools.BedTool(d_files['genome_masks'])

    tree = simulation.filter_masked_sites(tree, sim_data, gmask)

    d_fst = fst(tree)

    d_fst['sim_id'] = simid
    d_fst['category'] = vcat

    return d_fst


def main(tree_files, ncores, outfile):
    """
    main function 
    Args:
        tree_files a list of tskit tree files
        ncores: No. of cores for parallel processing
        outfile: table to save results
    """
    with mp.Pool(ncores) as pool:
        results = pool.map(get_fsts, tree_files)

    pd.concat(results).to_csv(outfile, index=False)


if __name__ == "__main__":
    *arg_tree_files, arg_ncores, arg_outfile = sys.argv[1:]
    arg_ncores = int(arg_ncores)
    main(arg_tree_files, arg_ncores, arg_outfile)
