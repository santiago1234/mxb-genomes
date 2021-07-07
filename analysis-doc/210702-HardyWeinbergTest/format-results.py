"""
Make a nice table with the results.
"""

import pandas as pd
import numpy as np
import glob
import os

spectrums = glob.glob("results/sfs*.csv")


def get_data_from_filename(sfn):
    """
    get the info from filename
    Note: the input sfn should look like: sfs-hw-e-6-nopasan-MX.csv
    """
    base = os.path.basename(sfn)
    d = base.split('-')
    hw_pass = d[2]
    sample = d[-1][:-4]
    return hw_pass, sample


def load_spectrum(sfn):
    spec = pd.read_csv(sfn)
    spec = (
        spec
        .rename({'index': 'pop1'}, axis=1)
        .melt(id_vars=['pop1'], value_name='freq', var_name=['pop2'])
    )
    spec[['pop1', 'n1']] = spec.pop1.str.split('-', 1, expand=True)
    spec[['pop2', 'n2']] = spec.pop2.str.split('-', 2, expand=True)

    # add the data
    spec[['hw_pass', 'Sample']] = get_data_from_filename(sfn)
    return spec


spectrums = pd.concat([load_spectrum(x) for x in spectrums])
spectrums.to_csv("results/all-spectrums-tidy.csv", index=False)


# the final spectrum

sfn = "results/final-sfs.csv"
spec = pd.read_csv(sfn)
spec = (
    spec
    .rename({'index': 'pop1'}, axis=1)
    .melt(id_vars=['pop1'], value_name='freq', var_name=['pop2'])
)
spec[['pop1', 'n1']] = spec.pop1.str.split('-', 1, expand=True)
spec[['pop2', 'n2']] = spec.pop2.str.split('-', 2, expand=True)
spec.to_csv("results/final-sfs-tidy.csv")
