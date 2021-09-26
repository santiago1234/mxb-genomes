"""
Compute Joint Site Frequency Spectrum for a set of biallelic variants.
author: santiago medina
help: python jsfs.py --h 
"""
import numpy as np
import pandas as pd
import argparse
from Bio import SeqIO
import moments
import allel
import pickle


# vcf = "../results/data/210713-HardyW-filters/1TGP_and_50MXB-chr22-snps-vep-mask-HW-GRCh38.vcf.gz"
# poplabs = "poplabs.csv"
# ancgenome = '/Users/santiagomedina/tmp/ancestral-genome-autosomes.fasta'
# outprefix = "tmp"


def load_data(vcf, ancgenome, poplabs):
    """
    Args:
    vcf:
        vcf: str path to vcf file
        anc: str path to fasta file .fa or .fasta
        poplabs: str path to csv file mapping pop labels
            to samples
    Returns:
        vcf: 
    """
    vcf = allel.read_vcf(
        vcf, alt_number=1)  # alt_number set to 1 (only biallelic SNPs)
    poplabs = pd.read_csv(poplabs)
    # check correct format of poplabs table
    if 'Samplename' not in poplabs.columns:
        raise ValueError('Samplename column not in poplabs')

    if 'Population' not in poplabs.columns:
        raise ValueError('Population column not in poplabs')

    # check all samples in poplabs present in vcf
    for s in poplabs['Samplename'].tolist():
        if s not in vcf['samples']:
            raise ValueError('sample {} not found in vcf file'.format(s))
    # The function Bio.SeqIO.to_dict() will use the
    # record ID as the dictionary key
    ancgenome = SeqIO.to_dict(SeqIO.parse(ancgenome, "fasta"))

    # check chromosomes in VCF are in ancestral genome
    chrs_vcf = [x for x in np.unique(vcf['variants/CHROM'])]

    for c in chrs_vcf:
        if c not in ancgenome.keys():
            raise ValueError(
                'chromosome {} not found in ancestral genome'.format(c))

    return vcf, poplabs, ancgenome


def count_allels(vcf, poplabs):
    """
    Return:
        dict mapping populations to allel counts array
    """
    ga = allel.GenotypeArray(vcf['calldata/GT'])
    # make dic mapping populations to Samplenames
    subpops = (
        poplabs
        .groupby('Population')
        .Samplename
        .apply(list)
        .to_dict()
    )
    # we need to map the pop to indices
    samples = list(vcf['samples'])

    def indexes_in_vcf(sl):
        return [samples.index(s) for s in sl]

    subpops = {pop: indexes_in_vcf(subpops[pop]) for pop in subpops.keys()}
    return ga.count_alleles_subpops(subpops, max_allele=1)


def annotate_ancestral_allel(vcf, ancgenome):
    """
    annotate the ancestral allele
    """
    ancgenome_str = {chrn: str(ancgenome[chrn].seq)
                     for chrn in ancgenome.keys()}

    def get_anc_allele(chrn, position):
        # subtract -1 to match position in sequence
        return ancgenome_str[chrn][position - 1]

    aa = [get_anc_allele(c, p) for (c, p) in zip(
        vcf['variants/CHROM'], vcf['variants/POS'])]

    d = {
        'ID': vcf['variants/ID'],
        'REF': vcf['variants/REF'],
        'ALT': vcf['variants/ALT'],
        'AA': aa
    }
    return pd.DataFrame(d)


def ancestral_allel_category(aa_table):
    """
    Generate categories for the ancestral allele
    Returns:
        dict with information
    """
    stats = dict()
    stats['total'] = aa_table.shape[0]
    stats['lowercase'] = aa_table.AA.isin(list('acgt')).values
    stats['uppercase'] = aa_table.AA.isin(list('ACGT')).values
    stats['N'] = (aa_table.AA == 'N').values
    stats['insertion_in_extant'] = (aa_table.AA == '-').values
    stats['no_coverage'] = (aa_table.AA == '.').values
    stats['REF_is_AA'] = (aa_table.AA == aa_table.REF).values
    stats['ALT_is_AA'] = (aa_table.AA == aa_table.ALT).values
    not_ref_alt = np.logical_and(
        np.logical_not(stats['REF_is_AA']),
        np.logical_not(stats['ALT_is_AA'])
    )
    stats['not_REF_or_ALT'] = np.logical_and(
        stats['uppercase'], not_ref_alt)
    return stats


def summary_aastats(stats):
    """
    Returns a string summarizing the ancestral allele information
    """
    discarded = sum(
        [stats[x].sum() for x in ['lowercase', 'insertion_in_extant',
                                  'no_coverage', 'not_REF_or_ALT', 'N']]
    )
    kept = stats['REF_is_AA'].sum() + stats['ALT_is_AA'].sum()
    summary = '''Total SNPs: {0}


    SNPs discarded: {1}, {11} %
        Why?:
            - low confident (acgt): {2}
            - failure (N): {3}
            - extant species contains insertion (-): {4}
            - no coverage in aligment (.): {5}
            - ancestral high confidence but is not ALT or REF: {6}


    SNPs kept (only high confidence: ACGT): {7},  {8} %
            - REF is AA: {9}
            - ALT is AA: {10}

    '''.format(
        stats['total'],
        discarded,
        stats['lowercase'].sum(),
        stats['N'].sum(),
        stats['insertion_in_extant'].sum(),
        stats['no_coverage'].sum(),
        stats['not_REF_or_ALT'].sum(),
        kept,
        round(100 * kept / stats['total'], 3),
        stats['REF_is_AA'].sum(),
        stats['ALT_is_AA'].sum(),
        round(100 * discarded / stats['total'], 3)
    )
    print(summary)
    return summary


def locate_variant_in_sfs(allele_counts, pops_to_index, var_index):
    """
    Where the variant falls in the SFS?
    Bassed on the derived allel count.
    """

    var_index = int(var_index)
    n_pops = len(pops_to_index.keys())

    pops_in_index_order = [pops_to_index[d] for d in range(n_pops)]

    var_location = [allele_counts[p][var_index][1]
                    for p in pops_in_index_order]
    return tuple(var_location)


def initialize_jsfs(poplabs):
    """
    Initialize the joint SFS
    """
    haploid_size = poplabs.groupby('Population').Population.count().to_dict()
    ndim = len(haploid_size.keys())

    def axis_length(hapsize):
        return hapsize * 2 + 1

    pops_to_index = {k: i for (i, k) in enumerate(haploid_size.keys())}
    index_to_pops = {i: k for (i, k) in enumerate(haploid_size.keys())}

    dimesion = [axis_length(haploid_size[index_to_pops[x]])
                for x in range(ndim)]
    sfs = np.zeros(dimesion, dtype=np.int32)
    return sfs, index_to_pops, pops_to_index


def populate_sfs_with_variants(sfs, allele_counts,  index_to_pops):
    n_vars = allele_counts[index_to_pops[0]].shape[0]
    for K in range(n_vars):
        var_loc = locate_variant_in_sfs(allele_counts, index_to_pops, K)
        sfs[var_loc] += 1
    return sfs


def polarize_counts(counts):
    """
    Counts: AlleleCountsArray
    switch alternative counts for reference counts
    """
    polarized = counts.copy()
    polarized[:, 0] = counts[:, 1]
    polarized[:, 1] = counts[:, 0]
    return polarized


def run(args):
    vcf = args.vcf
    poplabs = args.poplabs
    ancgenome = args.ancgenome
    outprefix = args.out

    print('loading data ...')
    vcf, poplabs, ancgenome = load_data(vcf, ancgenome, poplabs)

    print('extracting ancestral alleles ...')
    aa_table = annotate_ancestral_allel(vcf, ancgenome)

    print('computing jSFS ...')
    allele_counts = count_allels(vcf, poplabs)
    aa_condition = ancestral_allel_category(aa_table)
    stats = summary_aastats(aa_condition)
    sfs, index_to_pops, pops_to_index = initialize_jsfs(poplabs)

    # Polarize ancestral allele
    # REF = AA
    ac_ref_aa = {
        pop: allele_counts[pop][aa_condition['REF_is_AA']]
        for pop in allele_counts.keys()
    }
    sfs = populate_sfs_with_variants(sfs, ac_ref_aa, index_to_pops)

    ac_alt_aa = {
        pop: allele_counts[pop][aa_condition['ALT_is_AA']]
        for pop in allele_counts.keys()
    }

    ac_alt_aa_polarizedd = {pop: polarize_counts(
        ac_alt_aa[pop]) for pop in allele_counts}

    sfs = populate_sfs_with_variants(sfs, ac_alt_aa_polarizedd, index_to_pops)

    # order pops by index
    pop_ids = [index_to_pops[i] for i in range(len(index_to_pops))]
    spectrum = moments.Spectrum(sfs, pop_ids=pop_ids, data_folded=False)

    # save output
    stats_file = open(outprefix + '-stats.txt', "w")
    stats_file.write(stats)
    stats_file.close()

    spec_file = open(outprefix + '-spectrum.pkl', 'wb')
    pickle.dump(spectrum, spec_file)
    spec_file.close()


def main():
    parser = argparse.ArgumentParser(
        description="Compute joint SFS from vcf for a set of populations")

    parser.add_argument(
        "-vcf",
        help="VCF file. May be uncompressed or gzip-compatible  compressed file.",
        dest="vcf", type=str, required=True)

    parser.add_argument(
        "-poplabs",
        help="""
        Populations Information.
        csv file with columns Samplename (for samples in vcf) and
        Population.
        If there are k populations in the column Population the SFS will be k-dimensional.
        """,
        dest="poplabs", type=str, required=True)

    parser.add_argument(
        "-ancgenome",
        help="""
        Fasta file with ancestral genome. IDs in fasta should match chromosomes
        in VCF file.
        """,
        dest="ancgenome", type=str, required=True)

    parser.add_argument(
        "-out",
        help="""
        outprefix for output files. Two output files are written
         {out}-stats.txt summarizing ancestral allele state and 
         {out}-spectrum.pkl which contains the SFS in a moments.Spectrum object.
        """,
        dest="out", type=str, default="sfs"
    )
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
