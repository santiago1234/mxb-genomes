import numpy as np
import pandas as pd
import allel

# load files
vcf_file = '../../../../../tmp/1TGP_and_50MXB-chr22-snps-GRCh38.vcf.gz'
pop_info_file = 'data/popinfo.txt'
spectrum_out_file = 'data/sfs_subpops.csv'

# only load used fields from vcf
vcf = allel.read_vcf(vcf_file, fields=['samples', 'calldata/GT'])
pops = pd.read_csv(pop_info_file)

# get list of population names
populations = np.unique(pops.Population.to_numpy())

def get_population_indices(vcf, population, pop_info):
    """
    This function returns the indices in the vcf samples
    that are from the sampe population
    Params:
        vcf, vcf with samples
        population, str the population name
        pop_info: frame with columns Sample & Population
    """
    samples_in_population = pops[pops.Population == population].Sample.to_list()
    indices = [i for i, value in enumerate(vcf['samples']) if value in samples_in_population]
    
    return indices

subpops_indices = {pop:get_population_indices(vcf, pop, pops) for pop in populations}

# Allele counts per populations
# I use max_allele to only consider biallelic snps
ga = allel.GenotypeArray(vcf['calldata/GT'])
ac_subpops = ga.count_alleles_subpops(subpops_indices, max_allele=1)

def get_sfs_folded(sub_pop):
    """
    obtains the folded site frequency spectrum
    Returns: pd.DataFrame with columns:
        Population: name for population
        k: minor allele count
        sfs: frequency
    """
    sfs_f_subp = allel.sfs_folded(ac_subpops[sub_pop])
    k = list(range(sfs_f_subp.size))
    sfs_data = {"Population": sub_pop, "k": k, 'sfs': sfs_f_subp}
    
    return pd.DataFrame(sfs_data)

sfss = [get_sfs_folded(s) for s in ac_subpops.keys()]
sfs_all_subpops = pd.concat(sfss)
sfs_all_subpops.to_csv(spectrum_out_file, index=False)
