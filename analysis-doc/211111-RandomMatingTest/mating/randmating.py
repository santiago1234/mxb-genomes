import pandas as pd
import random


def parental_passed_genome(parent_data, chromosomes):
    """
    Args:
        parent_data: tracts data with haplotype information. This
            is a table containing a column for chrn (chromosome) and
            haplo (haplotype).
        chromosomes: list of chromosome names.

    Each parent has two chromosome copies (A or B), for each
    chromosome that will be passed to the next generation randomly
    select one copy.
    """
    n = len(chromosomes)
    haps_inherited = random.choices(['A', 'B'], k=n)
    inherited_chrs = pd.DataFrame({'chrn': chromosomes, 'haplo': haps_inherited})

    inherited_chrs = pd.merge(inherited_chrs, parent_data)
    # elimina la columna haplo
    return inherited_chrs.drop(['haplo'], axis=1)


def newborn(tractspop, chromosomes=list(range(1, 23))):
    """
    Args:
        tractspop: dict mapping individuals to tracts dta.
    """
    # make a random id
    id_new = ''.join(random.choices([str(x) for x in range(9)], k=4))
    # select parents at random
    mom, dad = random.sample(tractspop.keys(), 2)
    id_new = 'id' + id_new + '_' + mom + '_X_' + dad

    # recombination and mating
    retrive_tracts = lambda x: pd.concat(tractspop[x])
    mom_genome = parental_passed_genome(retrive_tracts(mom), chromosomes)
    dad_genome = parental_passed_genome(retrive_tracts(dad), chromosomes)

    return id_new, mom_genome, dad_genome
