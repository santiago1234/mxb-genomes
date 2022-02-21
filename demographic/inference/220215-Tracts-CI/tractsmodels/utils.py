import os
import pandas as pd
from scipy.stats import poisson
import glob
import numpy as np

_defaul_labs = ['NAT', 'EUR', 'AFR']

LABELS = {
    'CLM': _defaul_labs,
    'MXL': _defaul_labs,
    'PEL': _defaul_labs,
    'PUR': ['AFR', 'EUR', 'NAT']
}


# string between individual label and haploid chromosome id in input file
inter = "_anc"
# string at the end of input file. Note that the file should end in ".bed"
end = "_cM.bed"


def list_individuals_in_dir(directory, inter='_anc', end='_cM.bed'):
    """
    Args:
        directory: Path to directory with tracts input bed files.
        inter: string between individual label and haploid chromosome id in input file
        end: string at the end of input file. Note that the file should end in ".bed"
    """
    _files = os.listdir(directory)
    files = [filex
             for filex in _files
             if filex.split('.')[-1] == "bed"]
    ind_names = list(set(file.split('_')[0] for file in files))
    return ind_names


def bootsamp(num):
    '''generates a list of positions of the samples to pick in a bootstrap'''
    return np.random.choice(range(num), replace=True, size=num)


def get_bootstrap_replicata(pop, seed):
    if seed == 0:
        return pop
    else:
        np.random.seed(seed)
        indivs = pop.indivs
        bootorder = bootsamp(len(indivs))
        indivs2 = [indivs[i] for i in bootorder]
        pop.indivs = indivs2
        return pop


# helper functions to load tracts data results

def _list_tracts_files(path_to_files, mdl, population, bootstrap):
    """
    list tracts output files
    """
    file_basename = f'{path_to_files}/{population}-{mdl}-boot{bootstrap}'
    tfs = glob.glob(file_basename + '*')
    if len(tfs) == 0:
        raise ValueError(
            f'Files not found. Pattern: {file_basename} did not match anything.')

    return tfs


def _tfile_getter(file_list, termination):
    """
    Returns the file that ends with termination
    """
    try:
        fp, *_ = [x for x in file_list if x.endswith(termination)]
        return fp
    except ValueError:
        print(f'not files with termination: {termination}')


def confIntMeanPoisson(mu, conf=0.8):
    return poisson.interval(alpha=conf, mu=mu)


def _read_tf(file_list, termination):
    tf = _tfile_getter(file_list, termination)
    return np.loadtxt(tf)



def load_boot_rest(tfs, ancestries):
    """
    Args:
        tfs: list of tracts files
        ancestries: list, the population labels, same order as
        it was used to fit the model.
    """

    bins = _read_tf(tfs, 'bins')
    dat = _read_tf(tfs, 'dat')
    pred = _read_tf(tfs, 'pred')

    def make_frame(arr_d, var_name):
        "helper function to construc data frame"
        df = pd.DataFrame(arr_d.T, columns=ancestries)
        df['bins'] = bins
        df = df.melt(id_vars=['bins'],
                     var_name="Ancestry", value_name=var_name)
        return df

    dat = make_frame(dat, 'dat')
    pred = make_frame(pred, 'pred')
    data = pd.merge(dat, pred, how='inner', on=['bins', 'Ancestry'])
    c_i = data.pred.map(confIntMeanPoisson)
    data['ci_l'] = c_i.map(lambda x: x[0])
    data['ci_u'] = c_i.map(lambda x: x[1])
    return data


def _get_ancestries(tfs):
    """
    tfs: list of tracts files
    """
    labels = _tfile_getter(tfs, 'tsv')
    labels = pd.read_csv(labels, sep='\t')
    ancestries = labels.Anc.tolist()
    return ancestries


def ancestry_data_with_fits(path_to_files, mdl, population, bootstrap):
    """
    Returns a data frame with the following columns:
        bins: bins (cM)
        Ancestry: The Ancestry labels.
    """
    # tracts files list
    tfs = _list_tracts_files(path_to_files, mdl, population, bootstrap)
    fit_data = load_boot_rest(tfs, _get_ancestries(tfs))
    # Add the metadata to the table
    fit_data['Population'] = population
    fit_data['mdl'] = mdl
    fit_data['bootstrap'] = bootstrap
    return fit_data

