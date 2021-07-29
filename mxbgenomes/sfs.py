import allel
import moments
import numpy as np
import pandas as pd


# test_data

vcf_file = "/Users/santiagomedina/tmp/1TGP_and_50MXB-chr22-snps-vep-GRCh38.vcf.gz"
aa_file = "/Users/santiagomedina/tmp/aa-chr22.csv"


def ancestal_allel_stats(df):
    """
    Compute the ancestral allel statistics
    Args:
        df: DataFrame where each row is a variant
        and has REF, ALT, and AA columns
    Returns:
        stats: stats, how many REF is ALT?, etc.
    """
    aa_missing = df.AA.isnull().to_numpy()
    aa_is_dot = (df.AA == '.').to_numpy()
    aa_is_N = (df.AA == "N").to_numpy()  # N in fasta file
    aa_is_ref = (df.AA == df.REF).to_numpy()
    aa_is_alt = (df.AA == df.ALT).to_numpy()
    aa_unknown = (aa_is_N + aa_is_dot + aa_missing)
    # the next case is when the ancestral allel exits
    # but it is neither the alternative or the reference
    aa_not_ref_or_alt = np.logical_and(np.logical_and(~aa_is_ref, ~aa_is_alt), ~aa_unknown)
    stats = {
        'case': ['aa_is_ref', 'aa_is_alt', 'aa_is_dot', 'aa_is_N', 'aa_missing', 'aa_not_ref_or_alt'],
        'n': [aa_is_ref, aa_is_alt, aa_is_dot, aa_is_N, aa_missing, aa_not_ref_or_alt]
    }
    stats = pd.DataFrame(stats)
    stats['n'] = stats['n'].map(np.sum)
    stats['p'] = stats['n'] / stats['n'].sum() * 100
    return stats


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
        stats: pd.DataFrame, statistics for the ancestral allel and the vcf file
    """
    vcf = allel.read_vcf(
        vcf_file, alt_number=1)  # set alt_number=1 for biallelic variants
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
    # stats for th ancestral allel
    stats = ancestal_allel_stats(df)

    # Sanity check, that the variants are aligned with the vcf

    if not np.all(vcf['variants/ID'] == df.index):
        raise ValueError('sanity check not passed, data is no aligned')
    # We want to switch the reference allel
    # for the alternative allel when the
    # alternative allel is the ancestral.
    is_alt_aa = (df.ALT == df.AA).to_numpy()

    # Extract the genotype data

    ga = vcf['calldata/GT']

    return ga, is_alt_aa, vcf, stats


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
    Compute unfolded site frequency spectrum. In case
    the supplied ancestral allele matches the alternative allele
    it switches reference for alternative.
    Args:
        vcf_file: str, path to vcf file. Assume the SNPs in the
            vcf file are biallelic
        aa_file: str, path to csv file with variant ID and ancestral allel
            See:
                https://github.com/santiago1234/mxb-genomes/tree/main/analysis-doc/210506-AncestralAlleleData
        subpops: dict, maps supopulation names to sample names.
            sample names should be present in the vcf_file.
        project_haplod_size: int, haploid sample size to project
            SFS to. If None, the SFS is not projected.
    Returns:
        dict, mapping population names to SFS.
        stats: pd.DataFrame, statistics for the ancestral allel and the vcf file
    """
    print('loading vcf ...')
    ga, is_alt_aa, vcf, stats = load_GA_and_aa(vcf_file, aa_file)
    print('fixing ancestral allel ...')
    ga_fixed = fix_ancestral_allel(ga, is_alt_aa)
    ga_fixed = allel.GenotypeArray(ga_fixed)

    print('computing SFS ...')
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
        print('projecting SFS ...')
        sfs = {x: proyect_sfs(sfs[x], project_haplod_size, x) for x in sfs.keys()}

    return sfs, stats


def sfs_to_frame(sfs):
    """
    Put the project sfs in a nice format tidy frame
    Args:
        sfs: dict, mapping populations to sfs. All the sfs(s)
        have the same dimensions
    """
    df = (
        pd.DataFrame(sfs)
        .reset_index()
        .rename(columns={'index': 'n'})
        .melt(id_vars=['n'], var_name='Population', value_name='Freq')
    )
    return df


def joint_sfs(vcf_file, aa_file, populations, popinfo):
    """
    Compute the Joint Site Frequency Spectrum. It also polarizes
    to the ancestral allele.
    Args:
        vcf_file: str, path to VCF file.
        aa_file: str, path to csv file with variant ID and ancestral allel
            See:
                https://github.com/santiago1234/mxb-genomes/tree/main/analysis-doc/210506-AncestralAlleleData
        populations: list of populations to compute joint spectrum on.
        popinfo: population info. A data frame that has columns:
            Samplename -> The sample names that should be present in the vcf file.
            Subpopulation -> The subpopulations the populations should be found here. It
            may also contain other populations.
    Returns:
        tupple: spectrum, pop_index,
            spectrum -> multidimensional numpy array that represents the joint SFS.
            pop_index -> dict mapping indexes in multidimensional array to populations.
    """
    print('loading vcf ...')
    ga, is_alt_aa, vcf, stats = load_GA_and_aa(vcf_file, aa_file)
    print('polarizing to ancestral allel ...')
    ga_fixed = fix_ancestral_allel(ga, is_alt_aa)
    ga_fixed = allel.GenotypeArray(ga_fixed)


    # # Generate an array mapping pop names to pop indices in vcf

    subpops = (
        popinfo
        [popinfo.Subpopulation.isin(populations)]
        .groupby('Subpopulation')
        .apply(lambda x: x.Samplename.to_list())
        .to_dict()
    )


    # make sure all samples are in vcf
    [_val_samples(x, vcf) for x in subpops.values()]

    # Now get the indices in vcf file for the populations
    subpops_indices = {pop: get_population_indices(
        vcf, subpops[pop]) for pop in subpops}

    # allele counts
    ac = ga_fixed.count_alleles_subpops(subpops_indices, max_allele=1)

    # index to populations, In the joint SFS which pop is index 0, 1, 2, etc.?
    index_to_pops = {i: pop for i, pop in enumerate(subpops.keys())}

    diploid_size = lambda pop: 2 * len(subpops[pop]) + 1

    dimensions = [diploid_size(index_to_pops[i]) for i in index_to_pops]
    n_variants, _ = ac[list(subpops.keys())[0]].shape
    jsfs = np.zeros(dimensions, dtype=np.int32)


    def variant_position(var_index):
        """
        Where does the variant falls in the jsfs
        TODO:
        """
        # get the allel counts for the given variants
        # i use the index to get the data from the allel counts array
        variant_allele_counts = {pop:ac[pop][var_index] for pop in subpops.keys()}
        # get now the alternative allele count
        alternative_count = {pop:variant_allele_counts[pop][-1] for pop in subpops.keys()}
        # now get the index position for this variant in jsfs
        return tuple(alternative_count[index_to_pops[i]] for i in index_to_pops.keys())


    print('computing jsfs ...')
    for i in range(n_variants):
        jsfs[variant_position(i)] += 1

    return jsfs, index_to_pops

