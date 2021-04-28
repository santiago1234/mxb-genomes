"""
Combine all the results into a single table
"""
import pandas as pd
import glob
from os import path


files = glob.glob("data/*-mask-chr*.csv")
outfile = 'data/densities-all.csv'

def load_file(fname):
    # was vcf generated before of after mask filter?
    bsn = path.basename(fname)
    mask = bsn.split('-')[0]
    dat = pd.read_csv(fname)
    dat['mask_filter'] = mask
    return dat

densities = [load_file(x) for x in files]
densities = pd.concat(densities)
densities.to_csv(outfile, index=False)
