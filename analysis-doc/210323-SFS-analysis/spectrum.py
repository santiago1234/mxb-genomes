import numpy as np
import pandas as pd
import allel
import moments



def get_population_indices(vcf, population, pop_info):
    """
    This function returns the indices in the vcf samples
    that are from the sampe population
    Args:
        vcf, vcf with samples
        population, str the population name
        pop_info: frame with columns Sample & Population
    Returns:
        list containing the indices in the vcf samples
        that are from the given population
    """
    samples_in_population = pops[pops.Population == population].Sample.to_list()
    indices = [i for i, value in enumerate(vcf['samples']) if value in samples_in_population]
    
    return indices


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
