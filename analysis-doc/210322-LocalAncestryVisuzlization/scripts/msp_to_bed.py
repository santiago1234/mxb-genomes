#! /usr/bin/env python
from sys import path
import argparse
path.append('../../')
from mxbgenomes.localancestry import processxgmixout


def run(args):
    # set the arguments
    msp_file = args.msp_file
    individual = args.individual
    outfile = args.outfile
    d_ind = processxgmixout.get_data_for_indvididual(msp_file, individual)
    d_ind.to_csv(outfile, index=False, sep="\t")


def main():
    parser = argparse.ArgumentParser(
        description="Convert msp to bed like for a given individual"
    )
    parser.add_argument(
        "-msp", help="path to msp file", type=str, required=True,
        dest="msp_file"
    )
    parser.add_argument(
        "-ind", help="Sample identifier for individual",
        dest="individual",
        type=str, required=True)
    parser.add_argument(
        "-out", help="output name for bed file with ancstry assigments",
        type=str, required=True, dest="outfile"
    )
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
