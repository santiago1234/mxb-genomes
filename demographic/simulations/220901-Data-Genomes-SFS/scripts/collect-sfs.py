"""
usage:
    python collect-sfs.py sfs1 sfs2 ... outfile
"""
import pandas as pd
import glob
import os
import pickle
import moments
import sys

sfs = sys.argv[1:-1]
outfile = sys.argv[-1]


def load_sf(sfs_f):
    """
    Args:
        sf_f: file to site frequency spectrum
    """

    with open(sfs_f, "rb") as f:
        sf = pickle.load(f)

    # get metadata and pop
    cat = os.path.basename(sfs_f).split('cat_')[1].split('-')[0]
    pop = os.path.basename(sfs_f).split('pop_')[1].split('.pkl')[0]

    if cat == 'lof':
        cat = cat.upper()

    sf = sf.fold()

    d = {
        'Freq': sf.data,
        'Mask': sf.mask,
        'minor_allel_freq': range(sf.size)
    }

    sf = pd.DataFrame(d)

    sf = sf[~sf.Mask]
    sf['Pop'] = pop
    sf['Mut_Type'] = cat

    return sf.loc[:, ['Pop', 'Mut_Type', 'minor_allel_freq', 'Freq']]


all_sfs = pd.concat([load_sf(x) for x in sfs])

all_sfs.to_csv(outfile, index=False)
