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


def load_data(vcf, ancgenome, poplabs):
    """
    Args:
    vcf:
        vcf: str path to vcf file
        anc: str path to fasta file .fa or .fasta
        poplabs: str path to csv file mapping pop labels
            to samples
    Returns:
        vcf, poplabs
    """
    vcf = allel.read_vcf(
        vcf, alt_number=1)  # alt_number set to 1 (only biallelic SNPs)
    poplabs = pd.read_csv(poplabs)
    # check all samples in poplabs present in vcf
    for s in poplabs.Samplename.tolist():
        if s not in vcf['samples']:
            raise ValueError('sample {} not found in vcf file'.format(s))
    return vcf, poplabs


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



allele_counts = count_allels(vcf, poplabs)
sfs, pops_to_index, index_to_pops = initialize_jsfs(poplabs)

K = 1e3


def locate_variant_in_sfs(allele_counts, pops_to_index, var_index):

    var_index = int(var_index)
    n_pops = len(pops_to_index.keys())

    pops_in_index_order = [pops_to_index[d] for d in range(n_pops)]

    var_location = [allele_counts[p][var_index][1] for p in pops_in_index_order]
    return tuple(var_location)



# NOTE: I need to increase the legnth of each axis by 1.
def populate_sfs_with_variants(sfs, allele_counts,  pops_to_index):
    n_vars = allele_counts[pops_to_index[0]].shape[0]
    for K in range(n_vars):
        var_loc = locate_variant_in_sfs(allele_counts, pops_to_index, K)
        sfs[var_loc] += 1
    return sfs



