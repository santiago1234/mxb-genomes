"""
Usage: python get_r2.py {vcf_file} {bed_file} {rec_map} {output_file}

Args:
    vcf_file: path to VCF file
    bed_file: path to BED file
    rec_map: path to recombination map file should corresponed to the VCF file
    output_file: path to output file, to save r2 stats.

Computes the sums of LD statistics for data from the given population (assumes
all samples in the VCF are from the same population as provided).
"""

import sys, os, pickle

import numpy as np
import moments.LD

vcf_file, bed_file, rec_map, output_file = sys.argv[1:]

# the chromosome is found in the region bed file
with open(bed_file) as f:
    chrom = f.readline().split()[0]

r_bins = [
    0, 1e-7, 1e-6, 2e-6, 3e-6, 5e-6,
    1e-5, 2e-5, 3e-5, 5e-5,
    1e-4, 2e-4, 3e-4, 5e-4,
    1e-3, 2e-3, 3e-3, 5e-3, 
    1e-2] # out to 1 cM (1e-2)

# we'll use some functions in LD parsing to get the genotype matrix, etc
positions, genotypes, counts, sample_ids = moments.LD.Parsing.get_genotypes(
    vcf_file,
    bed_file=bed_file,
    chromosome=chrom,
    use_h5=False,
)

G = genotypes.to_n_alt()

pos_rs = moments.LD.Parsing._assign_recombination_rates(
    positions, rec_map, map_name="cM", cM=True, report=True
)

# use Huff/Rogers approach for computing r^2
def rogers_huff_r2(g1, g2):
    """
    g1 and g2 are the vectors of genotypes (0, 1, 2)
    """
    C = np.cov(g1, g2)[0, 1]
    vA = np.var(g1)
    vB = np.var(g2)
    r = C / np.sqrt(vA * vB)
    return r ** 2

# we'll store all r^2 values for each bin, then take the average at the end
# by averaging over all regions
r2s = [[] for r in r_bins[:-1]]
for ii, (g1, p1) in enumerate(zip(G, pos_rs)):
    for g2, p2 in zip(G[ii + 1:], pos_rs[ii + 1:]):
        dist = p2 - p1
        if dist >= r_bins[-1]:
            break
        r2 = rogers_huff_r2(g1, g2)
        bin_idx = np.where(dist >= r_bins)[0][-1]
        r2s[bin_idx].append(r2)

# to take averages later, and to reduce storage size, we'll take the sum over all
# values in each bin, and also record the total number of observations in each bin.
# then later, the average values can be computed by taking the grand sum divided
# by the grand total ("grand" meaning added over all regions)

sums = [sum(_) for _ in r2s]
totals = [len(_) for _ in r2s]

data = {
    "bins": r_bins,
    "sums": sums,
    "tots": totals
}

with open(output_file, "wb+") as fout:
    pickle.dump(data, fout)
