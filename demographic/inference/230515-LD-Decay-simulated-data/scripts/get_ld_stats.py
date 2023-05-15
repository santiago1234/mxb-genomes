"""
Usage: python ld_stats_region.py {region} {population}

Computes the sums of LD statistics for data from the given population (assumes
all samples in the VCF are from the same population as provided).
"""

import moments.LD
import sys, os, pickle

region = sys.argv[1]
pop = sys.argv[2]

bed_file = f"data/intervals/sim_{region}.tsv"
print(f'running region {region} for population {pop}')


# the chromosome is found in the region bed file
with open(bed_file) as f:
    chrom = f.readline().split()[0]

vcf_file = f"data/vcfs/sim{region}_{pop}.vcf" # format vcf, with chrom and pop
rec_map = f"data/recomb_map/sim_{region}-rmap.tsv" # format chromosome

out_file = f"results/ld_stats/{pop}-region{region}-ld_stats.pkl"

r_bins = [
    0, 1e-7, 1e-6, 2e-6, 3e-6, 5e-6,
    1e-5, 2e-5, 3e-5, 5e-5,
    1e-4, 2e-4, 3e-4, 5e-4,
    1e-3, 2e-3, 3e-3, 5e-3, 
    1e-2] # out to 1 cM (1e-2)

ld_stats = moments.LD.Parsing.compute_ld_statistics(
    vcf_file,
    bed_file=bed_file,
    chromosome=chrom,
    rec_map_file=rec_map,
    map_name="cM",
    cM=True,
    r_bins=r_bins,
    use_h5=False,
    report=False,
)

with open(out_file, "wb+") as fout:
    pickle.dump(ld_stats, fout)
