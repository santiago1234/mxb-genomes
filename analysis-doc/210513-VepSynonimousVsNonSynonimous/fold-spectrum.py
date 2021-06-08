"""
Fold the SFS for each variant category.
"""
import pandas as pd
import numpy as np
from glob import glob

spectrum = pd.read_csv("results/sfs-impact_HIGH.csv")


def fold_spectrum(s_unfolded):
    """
    Fold the SFS
    """
    n = np.floor(s_unfolded.size / 2)
    n = int(n)
    if n % 2 == 0:
        s_rev = s_unfolded[n:]
        s_rev = s_rev[::-1]
        folded = s_unfolded[:(n + 1)] + s_rev
        folded[-1] = folded[-1] / 2
    else:
        s_rev = s_unfolded[n:]
        s_rev = s_rev[::-1]
        folded = s_unfolded[:n] + s_rev
    return folded


def fold_pop(pop_spectrum):
    pop_spectrum = pop_spectrum.sort_values(by="n", axis=0)
    folded = fold_spectrum(pop_spectrum.Freq.to_numpy())
    n = range(folded.size)
    pop = pop_spectrum.Population.to_list()[0]
    sfolded = {
            'n': n,
            'Population': pop,
            'Minor_alle_freq': folded
            }
    return pd.DataFrame(sfolded)


def fold_var_cat(spectrum_path):
    """
    spectrum_path: path to specturm, example: results/sfs-impact_HIGH.csv
    """
    spectrum = pd.read_csv(spectrum_path)
    outname = spectrum_path.split('.csv')[0] + '-folded.csv'

    folded = (
        spectrum
        .groupby('Population')
        .apply(fold_pop)
        .reset_index(drop=True)
    )
    folded.to_csv(outname, index=False)


# get the file names to fold the project the spectrum
sfs_files = [x for x in glob('results/sfs-*.csv') if ('folded' not in x)]
# fold
[fold_var_cat(x) for x in sfs_files]

