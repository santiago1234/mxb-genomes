#! /usr/bin/env python
import pandas as pd
import argparse
from os import path
import sys
sys.path.append('../../')
from mxbgenomes.localancestry.tractslength import computer_tract_len_dist


def sample_name_from_filename(bed):
    """
    gets the sample name (Individual) name from
    input file path
    Args:
        bed: str, file path
    """
    sample_name = (
        path.basename(bed)
        .replace('allchrn-', '')
        .replace('.tsv', '')
    )
    return sample_name


def run(args):
    # input params
    bed = args.bed
    individual = sample_name_from_filename(bed)
    outfile = args.outfile
    bed = pd.read_table(bed)
    grouping_vars = ['chrn', 'Haplotype']  # I keep Individual

    tracts = (
        bed
        .groupby(grouping_vars)
        .apply(computer_tract_len_dist)
        .reset_index(level=grouping_vars)
    )

    # Agregate the Frequency by chromosomes and haplotypes
    tracts = (
        tracts
        .groupby(['Ancestry', 'tract_length'])
        [['Frequency']]
        .sum()
        .reset_index()
    )

    tracts['Individual'] = individual
    tracts.to_csv(outfile, index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Compute ancestry tract length distribution"
    )
    parser.add_argument(
        "-bed", type=str, required=True, dest="bed",
        help="Input bed file with ancestry assigment, all chromosomes, for a sample"
    )
    parser.add_argument(
        "-out", type=str, required=True, dest="outfile",
        help="Output file name for tract length distribution (csv)"
    )
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
