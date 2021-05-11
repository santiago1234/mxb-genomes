"""
Summarize tract length distribution across individuals.
This script computes the mean frequency across individuals
by population.
    The populations that I will use for summary
    are:
        -PEL
        -PUR
        -CLM
        -MXL
"""

from fnmatch import fnmatch
import numpy as np
import argparse
from scipy.stats import sem
import pandas as pd
import os
import sys
sys.path.append('../../')
from mxbgenomes.utils import load_1tgp_metada

# inpute parmeters
def run(args):
    # arguments

    path_to_panel = args.path_to_panel
    path_to_tracts = args.path_to_tracts
    patern_tracts = "*-tracts.csv"
    outputfile = args.outputfile

    # Load 1TGP pop info
    panel = load_1tgp_metada(path_to_panel)
    panel.drop(columns=['Superpopulation', 'Gender'], inplace=True)

    # load tracts files
    tracts_files = [os.path.join(path_to_tracts, x)
                    for x in os.listdir(path_to_tracts) if fnmatch(x, patern_tracts)]
    tracts = map(pd.read_csv, tracts_files)
    tracts = pd.concat(tracts)

    tracts.rename({'Individual': 'Samplename'}, axis=1, inplace=True)
    tracts = tracts.merge(panel, how='inner', on='Samplename')

    grp_vars = ['Ancestry', 'tract_length', 'Subpopulation']

    tracts_sumary = (
        tracts
        .groupby(grp_vars)['Frequency']
        .agg([('freq_mean', np.mean), ('freq_sem', sem)])
        .reset_index()
    )
    tracts_sumary.to_csv(outputfile, index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Summarize Trackt length distribution across individuals"
    )
    parser.add_argument(
        "-tracts",
        help="the path to the dir containing the tract length distribution for each individual",
        type=str, required=True,
        dest="path_to_tracts"
    )
    parser.add_argument(
        "-panel",
        help="resources/1tgp-samples-meta-data/integrated_call_samples_v3.20130502.all.panel",
        dest="path_to_panel",
        type=str, required=True
    )
    parser.add_argument(
        "-out", help="output name for file with mean frequency tract length",
        type=str, required=True, dest="outputfile"
    )
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
