"""
Compute Joint Site Frequency Spectrum for a set of biallelic variants.
author: santiago medina
"""
import numpy as np
import pandas as pd
from Bio import SeqIO
import moments
import allel


vcf = "../results/data/210713-HardyW-filters/1TGP_and_50MXB-chr22-snps-vep-mask-HW-GRCh38.vcf.gz"
poplabs = "poplabs.csv"
ancgenome = '/Users/santiagomedina/tmp/ancestral-genome-autosomes.fasta'


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


def initialize_jsfs(poplabs):
    """
    Initialize the joint SFS
    """
    haploid_size = poplabs.groupby('Population').Population.count().to_dict()
    ndim = len(haploid_size.keys())

    pops_to_index = {k: i for (i, k) in enumerate(haploid_size.keys())}
    index_to_pops = {i: k for (i, k) in enumerate(haploid_size.keys())}

    dimesion = [2 * haploid_size[index_to_pops[x]] for x in range(ndim)]
    sfs = np.zeros(dimesion, dtype=np.int32)
    return sfs, index_to_pops, pops_to_index


def annotate_ancestral_allel(vcf, ancgenome):
    """
    annotate the ancestral allele
    """
    ancgenome_str = {chrn: str(ancgenome[chrn].seq)
                     for chrn in ancgenome.keys()}

    def get_anc_allele(chrn, position):
        return ancgenome_str[chrn][position]

    aa = [get_anc_allele(c, p) for (c, p) in zip(
        vcf['variants/CHROM'], vcf['variants/POS'])]

    d = {
        'ID': vcf['variants/ID'],
        'REF': vcf['variants/REF'],
        'ALT': vcf['variants/ALT'],
        'AA': aa
    }
    return pd.DataFrame(d)


def ancestral_allel_stats(aa_table):
    """
    statistics for the ancestral allel
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
    not_ref_alt = np.logical_or(
        np.logical_not(stats['REF_is_AA']), 
        np.logical_not(stats['ALT_is_AA'])
    )
    stats['not_REF_or_ALT'] = np.logical_and(
        stats['uppercase'], not_ref_alt)
    return stats


vcf, poplabs, ancgenome = load_data(vcf, ancgenome, poplabs)
aa_table = annotate_ancestral_allel(vcf, ancgenome)

allele_counts = count_allels(vcf, poplabs)
sfs, pops_to_index, index_to_pops = initialize_jsfs(poplabs)

K = 1e3


def locate_variant_in_sfs(allele_counts, pops_to_index, var_index):

    var_index = int(var_index)
    n_pops = len(pops_to_index.keys())

    pops_in_index_order = [pops_to_index[d] for d in range(n_pops)]

    var_location = [allele_counts[p][var_index][1]
                    for p in pops_in_index_order]
    return tuple(var_location)


# NOTE: I need to increase the legnth of each axis by 1.
def populate_sfs_with_variants(sfs, allele_counts,  pops_to_index):
    n_vars = allele_counts[pops_to_index[0]].shape[0]
    for K in range(n_vars):
        var_loc = locate_variant_in_sfs(allele_counts, pops_to_index, K)
        sfs[var_loc] += 1
    return sfs
