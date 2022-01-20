'''
author: Santiago G. Medina-Muñoz
Compute joint SFS (Does not polarize ancestral allele).

Usage:
    python jsfs-nonPolarized.py <input.vcf> <poplabs.cvs> <output.pkl>

Args:
    * <input.vcf> VCF file. May be uncompressed or gzip-compatible  compressed file.

    * <poplabs.csv> Populations Information. A csv file with columns Samplename (for samples in vcf) and
         Population. If there are k populations in the column Population the SFS will be k-dimensional.

    * <output.pkl> Output which contains the SFS in a moments.Spectrum object.

NOTE: If you want a polarized SFS use jsfs.py
'''
import sys
import numpy as np
import pandas as pd
import moments
import allel
import pickle

# => => => => => =>  ||||||||| <= <= <= <= <= <= <= <= <= <= <= <= <= <= <=
# => => => => => =>  Input params <= <= <= <= <= <= <= <= <= <= <= <= <= <=
# => => => => => =>  ||||||||| <= <= <= <= <= <= <= <= <= <= <= <= <= <= <=

vcf, poplabs, out_spectrum = sys.argv[1:]

# => => => => => =>  ||||||||| <= <= <= <= <= <= <= <= <= <= <= <= <= <= <=
# => => => => => =>  FUNCTIONS <= <= <= <= <= <= <= <= <= <= <= <= <= <= <=
# => => => => => =>  ||||||||| <= <= <= <= <= <= <= <= <= <= <= <= <= <= <=


def load_data(vcf, poplabs):
    """
    Args:
    vcf:
        vcf: str path to vcf file
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


# => => => => => =>  ||||||||| <= <= <= <= <= <= <= <= <= <= <= <= <= <= <=
# => => => => => =>  MAIN CODE <= <= <= <= <= <= <= <= <= <= <= <= <= <= <=
# => => => => => =>  ||||||||| <= <= <= <= <= <= <= <= <= <= <= <= <= <= <=

print('loading data ...')
vcf, poplabs = load_data(vcf, poplabs)


print('computing jSFS ...')
allele_counts = count_allels(vcf, poplabs)
sfs, index_to_pops, pops_to_index = initialize_jsfs(poplabs)
sfs = populate_sfs_with_variants(sfs, allele_counts, index_to_pops)

# order pops by index
pop_ids = [index_to_pops[i] for i in range(len(index_to_pops))]

print('saving spectrum ...')
# save spectrum
spectrum = moments.Spectrum(sfs, pop_ids=pop_ids, data_folded=False)

spec_file = open(out_spectrum, 'wb')
pickle.dump(spectrum, spec_file)
spec_file.close()

print('done ...')
