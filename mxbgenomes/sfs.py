import allel
import moments
import numpy as np
import pandas as pd


# test_data

vcf_file = "/Users/santiagomedina/mxb-genomes/analysis-doc/210404-AncestralAlleleSFS/variant_effect_output.vcf.gz"
aa_file = "/Users/santiagomedina/mxb-genomes/analysis-doc/210404-AncestralAlleleSFS/ancestral.csv"


def load_GA_and_aa(vcf_file, aa_file):
    """
    Loads the Genotype Array from the VCF file.
    And obtain information about the ancestral allel.
    Args:
        vcf_file: str, path to vcf file. Assume the SNPs in the
            vcf file are biallelic
        aa_file: str, path to csv file with variant ID and ancestral allel
            See:
                https://github.com/santiago1234/mxb-genomes/tree/main/analysis-doc/210506-AncestralAlleleData
    Returns:
        genotype array: The VCF genotype array,
            3d array variants, samples, and chromosomes.
        is_alt_aa: logical vector where True indicates
            that the ALT allele is the AA. This array is
            aligned to ga, that is the i-th variant in ge
            is the i element in is_alt_aa
        vcf: vcf file object
    """
    vcf = allel.read_vcf(
        vcf_file, alt_number=1)  # set alt_number=1 to consider biallelic variants
    aa = pd.read_csv(aa_file)
    aa = aa.drop_duplicates(subset=['ID'])

    # Make a frame
    df = {
        'ID': vcf['variants/ID'],
        'REF': vcf['variants/REF'],
        'ALT': vcf['variants/ALT']
    }
    df = pd.DataFrame(df)

    # Add the ancestral allel
    df = pd.merge(df, aa, how='left', on='ID')
    df = df.set_index('ID')

    # Sanity check, that the variants are aligned with the vcf

    if not np.all(vcf['variants/ID'] == df.index):
        raise ValueError('sanity check not passed, data is no aligned')
    # We want to switch the reference allel
    # for the alternative allel when the
    # alternative allel is the ancestral.
    is_alt_aa = (df.ALT == df.AA).to_numpy()

    # Extract the genotype data

    ga = vcf['calldata/GT']

    return ga, is_alt_aa, vcf


def zeros_to_ones_anc_viceversa(x):
    if x == 0:
        return 1
    elif x == 1:
        return 0
    else:
        return x


def fix_ancestral_allel(ga, is_alt_aa):
    """
    Fix the ancestral allel to be the reference allel.
    Args:
        ga: genotype array
        is_alt_aa: logical vector where True indicates
            that the ALT allele is the AA. This array is
            aligned to ga, that is the i-th variant in ge
            is the i element in is_alt_aa
    Returns:
        Fixed genotype array. Fixed means that if ALT is
        the AA then we swicht a 0 by a 1 and viceversa.
    """
    ga_fixed = ga.copy()
    # I use this function becuase the genotype array may contain
    # other values besides 0s and 1s. For example -1  for missing
    # data and 2 for a multiallelic variant. I do not want to change
    # this values.
    fix_f = np.vectorize(zeros_to_ones_anc_viceversa)
    ga_fixed[is_alt_aa, :, :] = fix_f(ga[is_alt_aa, :, :])
    return ga_fixed


def get_population_indices(vcf, samples):
    """
    This function returns the indices in the vcf samples
    that are from the given population (samples list).
    Params:
        vcf, vcf with samples
        samples: list of sample names
    """
    vcf_samples = vcf['samples'].tolist()
    indices = [vcf_samples.index(s) for s in samples]
    return indices


def _val_sample(sample, samples_in_vcf):
    if sample not in samples_in_vcf:
        raise ValueError('Sample: {} not found in vcf samples'.format(sample))


def _val_samples(samples, vcf):
    """
    Make sure list of supplied samples are
    present in vcf file
    """
    [_val_sample(x, vcf['samples']) for x in samples]


def proyect_sfs(sfs, n, pop_id):
    """
    Args:
        sfs: np.array, The frequency spectrum data
        n: diploid population size to project to
        pod_id: str, population name
    Returns:
        np.array project spectrum
    """
    # make moments.Spectrun object
    spectrum = moments.Spectrum(sfs, pop_ids=[pop_id], data_folded=False)
    diploid_pop_size = 2 * n
    projected_sp = spectrum.project([diploid_pop_size])
    return np.array(projected_sp)
    

def sfs_unfolded(vcf_file, aa_file, subpops=None, project_haplod_size=None):
    """
    Args:
        subpops: dict, maps supopulation names to sample names.
            sample names should be present in the vcf_file.
        project_haplod_size: int, haploid sample size to project
            SFS to. If None, the SFS is not projected.
    Returns:
        dict, mapping population names to SFS.
    """
    ga, is_alt_aa, vcf = load_GA_and_aa(vcf_file, aa_file)
    print('fixing ancestral allel')
    ga_fixed = fix_ancestral_allel(ga, is_alt_aa)
    ga_fixed = allel.GenotypeArray(ga_fixed)

    print('computing SFS')
    if subpops is None:
        ac = ga_fixed.count_alleles(max_allele=1)
        sfs = allel.sfs(ac[:, 1])
        sfs = {'Global': sfs}  # Make the output a dict, so it is consistent

    else:
        # make sure all samples are in vcf
        [_val_samples(x, vcf) for x in subpops.values()]
        # Now get the indices in vcf file for the populations
        subpops_indices = {pop: get_population_indices(
            vcf, subpops[pop]) for pop in subpops}
        ac = ga_fixed.count_alleles_subpops(subpops_indices, max_allele=1)
        sfs = {x: allel.sfs(ac[x][:, 1]) for x in subpops}

    if project_haplod_size is not None:
        print('projecting SFS')
        sfs = {x: proyect_sfs(sfs[x], project_haplod_size, x) for x in subpops}

    return sfs

# helper code
popinfo = load_populations_info("../")

# Generate an array mapping pop names to pop indices in vcf
subpops = (
    popinfo
    .groupby('Subpopulation')
    .apply(lambda x: x.Samplename.to_list())
    .to_dict()
)
