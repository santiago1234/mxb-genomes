# Phasing MXB 50 genomes with 1TGP high coverage genomes
# The first part of the pipeline consist

# Note: for the rules in this file I use
# the prefix mp (marge phase) to avoid name
# conflicts with outher rules

from os import path

path_to_1TGP = config['path_to_1TGP']
oneTGP_pops = config['oneTGP_pops']


rule mp_subset_pop_list:
    input:
        pop_meta = "resources/1TGP-samples-meta-data/igsr-1000genomes.tsv"
    output:
        "results/data/210319-Phase-MXB-with-1TGP/merge-genomes/pops-to-subset.txt"
    shell:
        """
        cut -f1,4 {input} |\
            grep -E '{oneTGP_pops}' |\
            cut -f1 >{output}
        """

# chromosomes names must be ranamed in
# order to merge both datasets
rule mp_file_to_rename_chromosomes:
    output:
        "results/data/210319-Phase-MXB-with-1TGP/merge-genomes/rename-chromosomes.txt"
    shell:
        """
        touch {output}
        for i in {{1..22}}
        do
          echo "chr$i $i" >>{output}
        done
        """


rule mp_oneTGP_get_biallelic_and_subset_pops:
    # subset the population from 1TGP data
    # filter biallelic snps
    # also changes the name of the chromosome from chrN to N
    # to  match the chromosome names in the MXB genomes
    input:
        vcf = path.join(path_to_1TGP, "CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrn}.filtered.shapeit2-duohmm-phased.vcf.gz"),
        pops = "results/data/210319-Phase-MXB-with-1TGP/merge-genomes/pops-to-subset.txt",
        chr_names = "results/data/210319-Phase-MXB-with-1TGP/merge-genomes/rename-chromosomes.txt"
    output:
        temp("results/data/210319-Phase-MXB-with-1TGP/merge-genomes/1TGP-chr{chrn}-biallelic-subset.vcf.gz"),
        temp("results/data/210319-Phase-MXB-with-1TGP/merge-genomes/1TGP-chr{chrn}-biallelic-subset.vcf.gz.tbi")
    shell:
        """
        bcftools view -m2 -M2 -v snps {input.vcf} |\
            bcftools view -S {input.pops} |\
            bcftools annotate --rename-chrs {input.chr_names} -Oz -o {output[0]}
        bcftools index {output[0]} --tbi
        """


rule mp_50mxb_biallelic_snps:
    # SHAPEIT2 does not handle multiallelic variant phasing
    # this rule retrieves biallelic snp
    # and removes duplicated positions. For example chromosome 22 has
    # ~22 duplicated positions. This is a minimal number but shapeit is not
    # able to run with duplicated positions
    input:
        "results/data/raw-genomes/mxb-chr{chrn}-GRCh38.vcf.gz"
    output:
        temp("results/data/210319-Phase-MXB-with-1TGP/merge-genomes/mxb-unphased-chr{chrn}-GRCh38.vcf.gz"),
        temp("results/data/210319-Phase-MXB-with-1TGP/merge-genomes/mxb-unphased-chr{chrn}-GRCh38.vcf.gz.tbi")
    shell:
        """
        bcftools view -m2 -M2 -v snps {input} |\
            bcftools norm -d none -Oz -o {output[0]}
        bcftools index {output[0]} --tbi
        """

rule mp_merge_50G_with_1TG:
    # merge the 50 genomes with 
    # the 1TGP samples
    # I also set a variant id for variants
    # with missing id
    # note that I add in the bcftools pipeline
    # a new filter to keep biallelic variants
    # mergin the datasets my generate multiallelic variants
    input:
        vcf_50g = "results/data/210319-Phase-MXB-with-1TGP/merge-genomes/mxb-unphased-chr{chrn}-GRCh38.vcf.gz",
        index_50g = "results/data/210319-Phase-MXB-with-1TGP/merge-genomes/mxb-unphased-chr{chrn}-GRCh38.vcf.gz.tbi",
        vcf_1Tg = "results/data/210319-Phase-MXB-with-1TGP/merge-genomes/1TGP-chr{chrn}-biallelic-subset.vcf.gz",
        index_1Tg = "results/data/210319-Phase-MXB-with-1TGP/merge-genomes/1TGP-chr{chrn}-biallelic-subset.vcf.gz.tbi"
    output:
        "results/data/210319-Phase-MXB-with-1TGP/merge-genomes/1TGP_and_50MXB-chr{chrn}-snps-GRCh38.vcf.gz"
    shell:
        """
        bcftools merge --missing-to-ref \
            {input.vcf_1Tg} {input.vcf_50g} |\
             bcftools annotate --set-id \
             +'%CHROM\:%POS\:%REF\:%FIRST_ALT' |\
             bcftools view -m2 -M2 -v snps \
             -Oz -o {output}
        """


rule mp_phase_genomes:
    input:
        vcf = "results/data/210319-Phase-MXB-with-1TGP/merge-genomes/1TGP_and_50MXB-chr{chrn}-snps-GRCh38.vcf.gz",
        genetic_map = "resources/genetic-maps/chr{chrn}.b38.gmap"
    output:
        "results/data/210319-Phase-MXB-with-1TGP/phased/mxb-1tgp-chr{chrn}-snps-phased.sample",
        "results/data/210319-Phase-MXB-with-1TGP/phased/mxb-1tgp-chr{chrn}-snps-phased.haps"
    log:
        "results/logs/210319-Phase-MXB-with-1TGP/chrn{chrn}"
    params:
        out_basename = "results/data/210319-Phase-MXB-with-1TGP/phased/mxb-1tgp-chr{chrn}-snps-phased"
    threads: 15
    shell:
        """
        shapeit --input-vcf {input.vcf} -O {params.out_basename} \
            -T {threads} -M {input.genetic_map} --output-log {log}
        """


rule mp_convert_haps_to_vcf:
    input:
        "results/data/210319-Phase-MXB-with-1TGP/phased/mxb-1tgp-chr{chrn}-snps-phased.sample",
        "results/data/210319-Phase-MXB-with-1TGP/phased/mxb-1tgp-chr{chrn}-snps-phased.haps"
    output:
        "results/data/210319-Phase-MXB-with-1TGP/phased/mxb-1tgp-chr{chrn}-snps-phased.vcf.gz",
        "results/data/210319-Phase-MXB-with-1TGP/phased/mxb-1tgp-chr{chrn}-snps-phased.vcf.gz.tbi",
    params:
        basename = "results/data/210319-Phase-MXB-with-1TGP/phased/mxb-1tgp-chr{chrn}-snps-phased",
        out_vcf = "results/data/210319-Phase-MXB-with-1TGP/phased/mxb-1tgp-chr{chrn}-snps-phased.vcf"
    shell:
        """
        shapeit -convert --input-haps {params.basename} \
            --output-vcf {params.out_vcf}
        bcftools view {params.out_vcf} -Oz -o {output[0]}
        bcftools index {output[0]} --tbi
        rm {params.out_vcf}
        """


