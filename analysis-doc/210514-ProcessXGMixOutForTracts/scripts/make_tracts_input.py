"""
Takes the input bed file for one sample and generates the input format to run tracts
usage:
    python scripts/make_tracts_input.py <sample.bed> <samplename> <outdir>
"""
import os
import sys
import pandas as pd


bed_input = sys.argv[1]
samplename = sys.argv[2]
dest_dir = sys.argv[3]

bed = pd.read_table(bed_input)

def tracts_data(bed, hap):
    """
    Gets the data for one Haplotype and returns the order in tracts format
    """
    # the columns order that tracts expects
    col_order_tracts = ['chrn', 'spos', 'epos', 'Ancestry', 'sgpos', 'egpos']
    tracts = (
        bed[bed.Haplotype == hap]
        .loc[:, col_order_tracts]

    )
    return tracts


def out_name(samplename, dest_dir, hap):
    filename = samplename + '_anc_' + hap + '_cM.bed'
    filename = os.path.join(dest_dir, filename)
    return filename


tinput_A = tracts_data(bed, "A")
outname_A = out_name(samplename, dest_dir, "A")
tinput_A.to_csv(outname_A, index=False, sep="\t", header=False)


tinput_B = tracts_data(bed, "B")
outname_B = out_name(samplename, dest_dir, "B")
tinput_B.to_csv(outname_B, index=False, sep="\t", header=False)
