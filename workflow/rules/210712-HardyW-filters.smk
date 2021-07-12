"""
Pipeline to remove SNPs not in HW equilibrium.

======= >>>> 1st
First we remove SNPs not in HW equilibrium
within the continental populations:
    - AFR (YRI)
    - EUR (IBS + GBR)
    - EAS (CHB)
    - MXB (50 MXB)

We do not consider the admixed mexican population (MXL). We
use a relaxed p-value cut-off, 1e4.

======= >>>> 2nd
A second run of the HW equilibrium is to homogenize the 50 MXB
genomes with the 1TGP data. We have observed that there are
variants not called in the MXB genomes that were called in the 1TGP
and vice versa. This is likely due to differences in the variant
calling pipelines or also that the 50 MXB genomes were called in the
GRCh37 build (I lifter over to GRCh38) and the 1TGP in GRCh38.

For this I run HW equilibrium, in the following set of samples:
    - MX (50 MXB + MXL)
    - NAT (50 MXB + 1TGP-NAT)
        1TGP-NAT, is the samples from 1TGP that are a proxy for NAT ancestry.
        see: resources/1TGP-samples-meta-data/native-american.txt

Here, I use a strict HW pvalue of 0.05
"""
import pandas as pd
import sys
sys.path.append("./")
import os
from mxbgenomes.utils import load_populations_info
from os.path import join

# ***********************************************************
# ****** GLOBAL VARS ********
# ***********************************************************
# Define the p-value cutoffs
HWp_CONTINENTAL = 0.0001  # 1st
HWp_MX_NAT = 0.05  # 2nd

# The dir to save pipeline output
DESTINATION_DIR = "results/data/210713-HardyW-filters/"
CONTINENTAL = ['MXB', 'AFR', 'EUR', 'EAS']
MX_NAT = ['MX', 'NAT']

def subpops_to_samples():
    pops = load_populations_info("./")

    def get_continental_samples(pop):
        cont = pops[pops.Superpopulation == pop]
        cont = cont.Samplename.to_list()
        return cont


    mappings = {x: get_continental_samples(x) for x in CONTINENTAL}

    mxb_samples = get_continental_samples('MXB')
    mxl_samples = pops[pops.Subpopulation == 'MXL'].Samplename.to_list()
    nat1tgp_samples = (
        pd.read_table(
            "./resources/1TGP-samples-meta-data/native-american.txt",
            names=['Samplename']
        )
        .Samplename
        .to_list()
    )

    mappings['NAT'] = mxb_samples + nat1tgp_samples
    mappings['MX'] = mxb_samples + mxl_samples
    return mappings

SUBPOPS = subpops_to_samples()

# functions to apply on wildcars

def hw_pval(wildcars):
    """
    get the hw-pvalue bassed on HW.
    """
    if wildcars.subpop in CONTINENTAL:
        return HWp_CONTINENTAL
    else:
        return HWp_MX_NAT

def samples_list(wildcars):
    """
    Get the sample for the subpopulation.
    The output is a string were samples are separated by commas:
        A,B,C,D
    Intented for using with bcftools view -s
    """
    return ','.join(SUBPOPS[wildcars.subpop])


# *********************************************************************
# ********** PIPELINE CODE
# *********************************************************************

rule hw_all:
    input:
        expand("results/data/210713-HardyW-filters/1TGP_and_50MXB-chr{chrn}-snps-vep-mask-HW-GRCh38.vcf.gz", chrn=CHROMS)


rule hw_subpop_vcf:
    """
    Subset the sample for the subpopulation.
    Also remove non variant SNPs and keep only biallelic SNPs.
    """
    input:
        "results/data/210305-merged-with-1TGP/strict-mask/1TGP_and_50MXB-chr{chrn}-snps-vep-mask-GRCh38.vcf.gz"
    output:
        temp("results/data/210713-HardyW-filters/vcfs/subpop-{subpop}-chr{chrn}.vcf.gz")
    params:
        sl = samples_list
    shell:
        """
        bcftools view -s {params.sl} {input} |\
            bcftools view -c1 |\
            bcftools view -m2 -M2 -v snps -Oz -o {output}
        """


rule hw_test:
    """
    Run the HW and obtain p-values for each SNP.
    """
    input:
        "results/data/210713-HardyW-filters/vcfs/subpop-{subpop}-chr{chrn}.vcf.gz"
    output:
        temp("results/data/210713-HardyW-filters/hw/subpop-{subpop}-chr{chrn}.csv")
    message: "running HW test on {input}"
    conda:
        "../envs/scipy.yaml"
    shell:
        """
        python workflow/scripts/hw-test.py {input} {output}
        """

rule hw_snps_to_rm:
    """
    Get the SNPs that will be removed.
    """
    input:
        "results/data/210713-HardyW-filters/hw/subpop-{subpop}-chr{chrn}.csv"
    output:
        "results/data/210713-HardyW-filters/snps-to-remove/subpop-{subpop}-chr{chrn}.txt"
    params:
        hw_pval = hw_pval  # I am calling the function to get the cuttoff bassed on the population
    shell:
        """
        # Simply filter SNPs from the table
        awk -F, '{{ if ($3 < {params.hw_pval}) print $1}}' {input} >{output}
        """


rule hw_filter_VCF_by_HW:
    input:
        vcf = "results/data/210305-merged-with-1TGP/strict-mask/1TGP_and_50MXB-chr{chrn}-snps-vep-mask-GRCh38.vcf.gz",
        index = "results/data/210305-merged-with-1TGP/strict-mask/1TGP_and_50MXB-chr{chrn}-snps-vep-mask-GRCh38.vcf.gz.tbi",
        snps_to_drop = expand("results/data/210713-HardyW-filters/snps-to-remove/subpop-{subpop}-chr{{chrn}}.txt", subpop=CONTINENTAL + MX_NAT)
    output:
        all_snp_to_rm = "results/data/210713-HardyW-filters/snps-to-remove/all-chr{chrn}.txt",
        vcf = "results/data/210713-HardyW-filters/1TGP_and_50MXB-chr{chrn}-snps-vep-mask-HW-GRCh38.vcf.gz",
        index = "results/data/210713-HardyW-filters/1TGP_and_50MXB-chr{chrn}-snps-vep-mask-HW-GRCh38.vcf.gz.tbi"
    shell:
        """
        cat {input.snps_to_drop} | sort| uniq >{output.all_snp_to_rm}
        bcftools view --exclude ID==@{output.all_snp_to_rm} {input.vcf} -Oz -o {output.vcf}
        bcftools index {output.vcf} --tbi
        """

