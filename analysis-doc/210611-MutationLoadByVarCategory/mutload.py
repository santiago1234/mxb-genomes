import pandas as pd
import numpy as np


d = pd.read_csv("../210513-VepSynonimousVsNonSynonimous/results/sfs-missense.csv")

# Get the SFS
sfs = d[d.Population == "MXB"].Freq.to_numpy()

def l_fitness_reduction(q, h=0.25, s=0.08):
    """
    Fitnness reduction for allele
    Args:
        q: is the derived frequency of the mutation in population
        h: the dominance coefficient
        s: is the selection coefficient against the mutation
    """
    l = s * (2*h*q + (1 - 2*h)*q**2)
    return l

# Compute the load
# allele frequencies

diploid_size = len(sfs) - 1
allele_frqs = np.arange(0, diploid_size + 1) / diploid_size

# the load in each category 

x = (1 - l_fitness_reduction(allele_frqs))


def log_compute_load_sfs(sfs):
    diploid_size = len(sfs) - 1
    allele_frqs = np.arange(0, diploid_size + 1) / diploid_size
    fr = l_fitness_reduction(allele_frqs)
    fr = np.log(1 - fr)
    return fr * sfs, l_fitness_reduction(allele_frqs)


def spectrum_load(d_pop):
    d_pop = d_pop.sort_values('n')
    sfs = d_pop.Freq.to_numpy()
    d_pop['s_log_load'], d_pop['fitness_reduction'] = log_compute_load_sfs(sfs)
    return d_pop


L = d.groupby('Population').apply(spectrum_load).reset_index(drop=True)
L.to_csv('results/load-missense.csv', index=False)
