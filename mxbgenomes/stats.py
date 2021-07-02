"""
Some usefule statistics like HW test.
"""
import allel
import numpy as np
from scipy.stats import chisquare


def hw_test(ga):
    """
    Hardy Weinberg chisquare test. It is assumed only
    for biallelic loci.
    There is 1 degree of freedom (degrees of freedom for test for
    Hardy–Weinberg proportions are # genotypes − # alleles).
    Args:
        ga: allel.GenotypeArray
    Returns:
        chisq: The chi-squared test statistic.
        p: The p-value of the test.
    """
    # Compute the allel frequencies
    ac = ga.count_alleles(max_allele=1)
    ac_p = ac.to_frequencies()
    # get the number of genotypes from the allel counts
    n = ac.sum(axis=1) / 2  # the number of genotypes

    # Observed genotypes
    # In the genotype array 0 is the reference allele,
    # 1 is the first alternate allele.
    # compute genotype, 0 is AA, 1 is Aa, 2 is aa.
    genotypes = ga.sum(axis=2)
    AA_obs = (genotypes == 0).sum(axis=1)
    Aa_obs = (genotypes == 1).sum(axis=1)
    aa_obs = (genotypes == 2).sum(axis=1)
    observed = np.stack((AA_obs, Aa_obs, aa_obs), axis=1)

    # compute expected under HW
    # get the p and qs
    p = ac_p[:, 0]
    q = ac_p[:, 1]
    AA_exp = p**2 * n
    Aa_exp = 2 * p * q * n
    aa_exp = q**2 * n
    expected = np.stack((AA_exp, Aa_exp, aa_exp), axis=1)

    chi, p = chisquare(observed, expected, axis=1, ddof=1)
    return chi, p
