"""
Compute the 2-d SFS
usage:
    python scripts/generate_spectrums.py <NCORES>
"""
import sys
import multiprocessing as mp

import numpy as np
import pandas as pd

sys.path.append('../')
sys.path.append('../220910-Fst-Comparison/scripts/')

from simutils import utils, simulation
import fst_sim


def neutral_joint_spectrums(simid: int, rootpath='../../../'):
    '''
    Get the neutral (noncoding) joint spectrums
    Args:
        rootpath: path to project root dir mxb-genomes
    '''

    simdata, ts_sim, _ = utils.load_simulation_output(rootpath, simid)

    print('subsampling and adding neutral variation')
    ts_sim = simulation.subsample_individuals_pop(ts_sim, 50)
    sim_ts, _ = simulation.simulate_neutral_variation(ts_sim, simdata)

    return fst_sim.two_d_sfs(sim_ts, fold=False)


def combine_spectrums(spectrums):
    '''
    Combine the spectrums (sum them)
    '''

    def sum_spec(poppair):
        return sum(x[poppair] for x in spectrums)

    poppairs = spectrums[0].keys()
    poppairs = list(poppairs)

    combined = {pair: sum_spec(pair) for pair in poppairs}

    return combined


def gather_spectrums(up_to_n, ncores=10):
    """
    Get all spectrums form id 1 up_to_n
    """
    with mp.Pool(ncores) as pool:
        specs = pool.map(neutral_joint_spectrums, range(1, up_to_n + 1))

    return specs


def format_spectrum(spectrum):
    """
    Put the Spectrums in a nice pandas frame
    Args:
        spectrum: moments.Spectrum
    """

    pop1, pop2 = spectrum.pop_ids
    data  = spectrum.data

    # this are the fixed values, which we do not care
    data[0, 0] = np.nan
    data[-1, -1] = np.nan

    data = pd.DataFrame(data)
    data['Pop1_derived_freq'] = range(data.shape[0])

    data = data.melt(id_vars=['Pop1_derived_freq'], var_name='Pop2_derived_freq', value_name='Freq')
    data['pop_pair'] = f'{pop1}-{pop2}'
    return data


def main():
    n_cores = sys.argv[1]
    n_cores = int(n_cores)
    n_spectrums = 350
    specs = gather_spectrums(n_spectrums, ncores=n_cores)
    specs = combine_spectrums(specs)
    specs = pd.concat([format_spectrum(x) for x in specs.values()])
    specs.to_csv('results/joint-spectrums.csv')


if __name__ == "__main__":
    main()
