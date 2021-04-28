#! /usr/bin/env python
import argparse
import allel
import numpy as np
import pandas as pd


def variant_sites(gt):
    """
    Params:
        gt: Genotpye Array
    Returns:
        lgl array, True indicate that the particular
        site is a variant site. I define this if the
        heterozygositiy of the variant is more than 0
    """
    ga = allel.GenotypeArray(gt)
    af = ga.count_alleles().to_frequencies()
    het = allel.heterozygosity_expected(af, ploidy=2)
    return het > 0


def variant_density(variant_regions, min_pos, max_pos, window_size=100000):
    """
    Compute the variant density
    """
    bins = np.arange(min_pos, max_pos, window_size)

    # use window midpoints as x coordinate
    x = (bins[1:] + bins[:-1])/2

    # compute variant density in each window
    h, _ = np.histogram(variant_regions, bins=bins)
    y = h / window_size
    return x, y


def compute_var_density(callset, samples):
    """
    samples: lgl array, with True corresponding to the
    samples to use
    """
    gt = callset['calldata/GT']
    gt = gt[:, samples, :]  # subset the samples
    # get the variant position for these particular samples
    variants = variant_sites(gt)
    positions = callset['variants/POS']
    variants_positions = positions[variants]
    sp = positions.min()
    ep = positions.max()
    coord, density = variant_density(variants_positions, sp, ep)
    d = {
        'position_bp': coord,
        'variant_density': density
    }
    return pd.DataFrame(d)


def variant_density_mxb_and_1tgp(callset):
    """
    Compute the variant frequecy in the 50 MXB
    and in the 1TGP samples
    """
    # get the samples list
    samples = callset['samples']
    mxb_samples = np.array([x.startswith('MXB_') for x in samples])
    oneT_samples = np.logical_not(mxb_samples)

    density_50mxb = compute_var_density(callset, mxb_samples)
    density_oneT = compute_var_density(callset, oneT_samples)
    # add source info to merge tables
    density_50mxb['DataSource'] = "50 MXB"
    density_oneT['DataSource'] = "1TGP"
    # Normalize density by sample size
    density_50mxb['var_den_norm_sample'] = density_50mxb.variant_density / \
        mxb_samples.sum()
    density_oneT['var_den_norm_sample'] = density_oneT.variant_density / \
        oneT_samples.sum()
    return pd.concat([density_50mxb, density_oneT])


def run(args):
    # set the arguments
    vcf_file = args.vcf
    chrn = args.chr
    outfile = args.outfile
    callset = allel.read_vcf(vcf_file)
    var_den = variant_density_mxb_and_1tgp(callset)
    var_den['chromosome'] = chrn
    var_den.to_csv(outfile, index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Compute variant density in MXB and 1TGP samples"
    )
    parser.add_argument(
        "-vcf", help="path to VCF file", type=str, required=True,
        dest="vcf"
    )
    parser.add_argument(
        "-chr", help="chromosome name for the VCF file",
        dest="chr",
        type=str, required=True)
    parser.add_argument(
        "-out", help="output name for csv output table with variant density",
        type=str, required=True, dest="outfile"
    )
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
