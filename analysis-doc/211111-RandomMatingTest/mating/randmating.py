import pandas as pd
import random

mxl = load_tracts_pop('../../210514-ProcessXGMixOutForTracts/data/3-pops/tracts/MXL/')

population = mxl.keys()
chromosomes = list(range(1, 23))
## create a newborn sample

## sample parents without replacement

mom, dad = random.sample(population, 2)

# for each chromosome we randomnly choose one haplotype.

def parental_inherited_chromo(parent_data, chromosomes):
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
    inherited_chrs.drop(['haplo'], axis=1)
    return inherited_chrs


def newborn(tractspop):
    """
    Args:
        tractspop: dict mapping individuals to tracts dta.
    """
    # make a random id
    id_new = ''.join(random.choices([str(x) for x in range(9)], k=4))
    # select parents at random




