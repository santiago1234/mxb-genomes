"""
Module to load tracts results in a nice format
for visualization.
And some other functions.
"""
import pandas as pd
import os
import glob
import numpy as np
from scipy.stats import poisson


def _list_boostrap_files(dir_to_res, boot_n=0):
    """
    List the tracts output files for a particular boostrap run
    Returns:
        list of file for the given boostrap
    """
    start_pat = 'boot{n}_*'.format(n=boot_n)
    start_pat = os.path.join(dir_to_res, start_pat)
    return glob.glob(start_pat)


def _load_boot_file(dir_to_res, boot_n=0, filetype="dat"):
    """
    Args:
        filetype: str, either dat, bins, pred, etc.
    Loads the bins as an array
    """
    boot_file = [x for x in _list_boostrap_files(
        dir_to_res, boot_n) if x.endswith(filetype)][0]
    return np.loadtxt(boot_file)


def confIntMeanPoisson(mu, conf=0.8):
    return poisson.interval(alpha=conf, mu=mu + + 1e-9)


def load_boot_rest(dir_to_res, boot_n, labels):
    """
    Args:
        dir_to_res: the path to the directory with tracts output files
        labels: list, the population labels, same order as
        it was used to fit the model.
        boot_n: int, the boostrap number
    """
    bins = _load_boot_file(dir_to_res, boot_n, "bins")
    dat = _load_boot_file(dir_to_res, boot_n, "dat")
    pred = _load_boot_file(dir_to_res, boot_n, "pred")
    # Make data frame for dat

    def make_frame(arr_d, var_name):
        "helper function to construc data frame"
        df = pd.DataFrame(arr_d.T, columns=labels)
        df['bins'] = bins
        df = df.melt(id_vars=['bins'],
                     var_name="Ancestry", value_name=var_name)
        return df

    dat = make_frame(dat, 'dat')
    pred = make_frame(pred, 'pred')
    data = pd.merge(dat, pred, how='inner', on=['bins', 'Ancestry'])
    data['boot'] = boot_n
    c_i = data.pred.map(confIntMeanPoisson)
    data['ci_l'] = c_i.map(lambda x: x[0])
    data['ci_u'] = c_i.map(lambda x: x[1])
    return data


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

