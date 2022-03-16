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
    'PUR': ['EUR', 'NAT', 'AFR']
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


def _list_tracts_files4pops(path_to_files, mdl, bootstrap):
    """
    list tracts output files
    """
    file_basename = f'{path_to_files}/{mdl}-boot{bootstrap}'
    tfs = glob.glob(file_basename + '*')
    if len(tfs) == 0:
        raise ValueError(
            f'Files not found. Pattern: {file_basename} did not match anything.')

    return tfs


def ancestry_data_with_fits_4pops(path_to_files, mdl, bootstrap):
    """
    Returns a data frame with the following columns:
        bins: bins (cM)
        Ancestry: The Ancestry labels.
    """
    # tracts files list
    tfs = _list_tracts_files4pops(path_to_files, mdl, bootstrap)
    fit_data = load_boot_rest(tfs, _get_ancestries(tfs))
    # Add the metadata to the table
    fit_data['mdl'] = mdl
    fit_data['bootstrap'] = bootstrap
    return fit_data


# The next functions are helpfull for having a graphical representation
# of the model. See figure 5 in Gravel 2012.

def _ancestry_prop(pA_current, replacement, pA_r):
    """
    Args:
        pA_current: np.array, the current ancestry proportions in the population.
        replacement: fraction of new migrants entering the population. I assume
            this fraction will replace the current population.
        pA_r: np.array, propotion of migrants from each ancestry entering
            into the population.
    """
    # this operation follows from the law of total probability
    pA_new = pA_current * (1 - replacement) + replacement * pA_r
    return pA_new


def ancestry_prop_over_time(migmat, poplabels):
    """
    Reconstruct the ancestry proportion over time
    from the migration matrix.
    Args:
        migmat: np.array migration matrix. This migration matix
            is part of the Tracts output.
        poplabels: The population labels, the order must
            match the columns in the migration matrix.
    Returns: pd.DataFrame with columns:
        - ga: Generations ago from the present.
        - ancestry:
        - prop
    """
    # The last row of the migration matrix gives the initial
    # ancestry proportions
    current_p = migmat[-1, :]
    generations = list(range(migmat.shape[0]))
    generations = generations[::-1]

    anc_pt = np.zeros_like(migmat)
    anc_pt[generations[0], :] = current_p

    for t in generations[1:]:
        migrants = migmat[t, :]

        # Does this generation have a pulse of migration?
        # I do this to get the fraction of each ancestry
        # entering into the population.
        if migrants.sum() > 0:
            migrants_fraction = migrants / migrants.sum()
        else:
            migrants_fraction = np.zeros_like(migrants)

        current_p = _ancestry_prop(current_p, migrants.sum(), migrants_fraction)
        anc_pt[t, :] = current_p

    anc_pt = pd.DataFrame(anc_pt, columns=poplabels)
    anc_pt['ga'] = generations[::-1]
    return anc_pt


def migration_pulses(migmat, poplabels):
    """
    Get a summary of the migration pulses.
    Args:
        migmat: np.array migration matrix. This migration matix
            is part of the Tracts output.
        poplabels: The population labels, the order must
            match the columns in the migration matrix.
    """
    migmat = pd.DataFrame(migmat, columns=poplabels)
    migmat['ga'] = range(migmat.shape[0])

    pulses = migmat[migmat[poplabels].sum(axis=1) > 0.001]
    pulses['fr'] = migmat[poplabels].sum(axis=1)
    return pulses
