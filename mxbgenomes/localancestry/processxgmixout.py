"""
module to process XGMix output files
for downstream analysis
"""

import pandas as pd
import numpy as np
import os

_main_cols = ['#chm', 'spos', 'epos', 'sgpos', 'egpos', 'n snps']


def load_pop_codes(msp_file):
    """
    Reads the population codes from
    the msp_file. The population codes are included
    in the first line of the msp file.
    """
    mf = open(msp_file, 'r')
    header_codes = mf.readline()
    mf.close()
    pop_codes = (
        header_codes
        .replace('#Subpopulation order/codes: ', '')
        .strip()
        .split('\t')
    )
    return {int(b): a for (a, b) in [x.split('=') for x in pop_codes]}


def get_individual(msp, individual):
    """
    Retrieves the data for a particular
    individual
    Args:
        msp: pd.DataFrame
        individual: str, sample name
    """
    ind_haps = [individual + '.' + str(x) for x in [0, 1]]
    if ind_haps[0] not in msp.columns:
        raise ValueError('Supplied individual not found')
    return msp.loc[:, _main_cols + ind_haps]


def tidy_individual(ind_data):
    """
    Args:
        ind_data: pd.DataFrame
    """
    # make long format data
    ind_data_l = ind_data.melt(id_vars=_main_cols, value_name='Ancestry')
    # happlotype make 0s to A and 1s to B.
    # That way is more intuitive
    haplo_mapper = {'0': 'A', '1': 'B'}
    # the last digit is the haplotype
    ind_data_l['Haplotype'] = ind_data_l.variable.str[-1].map(haplo_mapper)
    ind_data_l['Individual'] = ind_data_l.variable.str[:-2]
    # rename columns to match Tracts readme columns
    col_renamer = {
        '#chm': 'chrn'
    }
    return ind_data_l.rename(col_renamer, axis=1)


def get_data_for_indvididual(msp_file, individual):
    """
    Args:
        msp_file: str, path to msp file generated by XGMIx
        individual: str, indentifier for individual
    """
    pop_codes = load_pop_codes(msp_file)
    msp = pd.read_table(msp_file, skiprows=1)
    ind_d = get_individual(msp, individual)
    ind_t = tidy_individual(ind_d)
    # change the pop codes to pop labels
    ind_t['Ancestry'] = ind_t.Ancestry.map(pop_codes)
    ind_t.drop('variable', axis=1, inplace=True)
    return ind_t


def list_msp_files(path_to_msp, fprefix="mdl-", fsufix=".msp.tsv"):
    """
    Get the list of the files generated by XGMix output.
    For all the autosomes.
    Args:
        path_to_msp: str, relative path to:
            /results/data/210409-local-ancestry/{3,4}-pops/predictions/
        Either for 3 or 4 populations.
    """
    files = [fprefix + str(x) + fsufix for x in range(1, 23)]
    files = [os.path.join(path_to_msp, x) for x in files]
    # make sure files exits
    for f in files:
        if not os.path.isfile(f):
            raise ValueError("file: {} not found".format(f))
    return files


def load_msp_file(msp_file):
    """
    Loads the msp_file as a pandas.DataFrame
    """
    msp = pd.read_table(msp_file, skiprows=1)
    p_codes = load_pop_codes(msp_file)
    samples_cols = [x for x in msp.columns if x not in _main_cols]
    msp[samples_cols] = msp[samples_cols].applymap(lambda x: p_codes[x])
    return msp


def load_xgmix_output(path_to_msp):
    """
    Load the XGMix ancestry assigments for
    all autosomes.
    Args:
        path_to_msp: str, relative path to:
            /results/data/210409-local-ancestry/{3,4}-pops/predictions/
        Either for 3 or 4 populations.
    Returns:
        pd.DataFrame
    """
    d = [load_msp_file(x) for x in list_msp_files(path_to_msp)]
    return pd.concat(d)



