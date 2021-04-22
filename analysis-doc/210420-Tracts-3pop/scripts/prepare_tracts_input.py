import pandas as pd
import os
import argparse
import sys
sys.path.append('../../')
from mxbgenomes.localancestry.tractslength import collapse_windows_to_tracts
# Format xgmix to run tracts


def get_haplo_data(tracts, hap):
    """
    Extract the data for one haplotype
    Args:
        tracts: pd.DataFrame,
        hap: str, either A or B for maternal or paternal haplotypes
    """
    cols_to_drop = ['len_bp', 'len_cm', 'Individual', 'Haplotype', 'n snps']
    # to keep the same order as the input files for tracts
    order_cols_tracs = ['chrn', 'spos', 'epos', 'Ancestry', 'sgpos', 'egpos']
    tracts_hap = (
        tracts[tracts.Haplotype == hap]
        .drop(labels=cols_to_drop, axis=1)
        .loc[:, order_cols_tracs]
    )
    return tracts_hap


def make_tracts_input(bed, sample_name, dest_dir):
    """
    Given the predicted Ancestry for a particular individual
    formats the file so that the output generated can be
    used to run Tracts analysis
    Args:
        bed: pd.DataFrame
        sample_name: str, sample name for Individual
        dest_dir: The dir to save output files: <Individual_anc_A_cm.bed>
            and <Individual_anc_B_cm.bed>
    """
    # Collapse to continuous Ancestry fragements
    tracts = (
        bed
        .groupby(['chrn', 'Haplotype'])
        .apply(collapse_windows_to_tracts)
        .reset_index(drop=True)
    )

    def make_filen(hap, sample_name):
        """Helper function to generate output file name"""
        filen = sample_name + '_anc_' + hap + '_cM.bed'
        return os.path.join(dest_dir, filen)

    d_hap_A = get_haplo_data(tracts, "A")
    d_hap_B = get_haplo_data(tracts, "B")
    d_hap_A.to_csv(make_filen("A", sample_name),
                   sep="\t", header=False, index=False)
    d_hap_B.to_csv(make_filen("B", sample_name),
                   sep="\t", header=False, index=False)


def run(args):
    bed = args.bed
    dest_dir = args.dest_dir
    bed = pd.read_table(bed)
    individual = bed.Individual.to_list()[0]
    make_tracts_input(bed, sample_name=individual, dest_dir=dest_dir)


def main():
    parser = argparse.ArgumentParser(
        description="Generate files with local ancestry for running Tracts"
    )
    parser.add_argument(
        "-bed",
        help="path to input bed file with ancestry information see: mxb-genomes/analysis-doc/210322-LocalAncestryVisuzlization/Snakefile",
        type=str,
        required=True,
        dest="bed"
    )
    parser.add_argument(
        "-ddir",
        help="destination dir, where to save output files?",
        dest="dest_dir",
        type=str, required=True
    )
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
