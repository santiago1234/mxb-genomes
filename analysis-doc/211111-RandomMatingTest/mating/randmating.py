import os
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
        tractspop: dict mapping individuals to tracts data.
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


def save_individual_tracts_to_file(id_new, mom_genome, dad_genome, outdir):
    """
    Save the data of an individual to a tracts (bed) format.
    """
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    file_A = os.path.join(outdir, id_new + '_anc_A_cM.bed')
    file_B = os.path.join(outdir, id_new + '_anc_B_cM.bed')

    mom_genome.to_csv(file_A, header=False, index=False, sep='\t')
    dad_genome.to_csv(file_B, header=False, index=False, sep='\t')
    pass


def simulate_random_maiting_pop(tractspop, outdir, n, seed):
    """
    Args:
        tractspop: dict mapping individuals to tracts data.
        seed: seed for random generation
    """
    random.seed(int(seed))
    [save_individual_tracts_to_file(*newborn(tractspop), outdir) for i in range(n)]
    pass



