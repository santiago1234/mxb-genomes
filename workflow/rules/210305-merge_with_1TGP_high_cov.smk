# Merge 50 genomes with 1TGP high coverage data
# note that I merge the phased bialleic snps from the 50 MXB genomes
# After I merge the data I will add the variant annotation with vep

from os import path

path_to_1TGP = config['path_to_1TGP']
oneTGP_pops = config['oneTGP_pops']

rule merge_all_50mxb_with_1TGP:
    input:
        expand("results/data/210305-merged-with-1TGP/1TGP_and_50MXB-chr{chrn}-snps-vep-GRCh38.vcf.gz", chrn=CHROMS)

rule subset_pop_list:
    #list with populations to subset from 1TGP
    input:
        pop_meta = "resources/1TGP-samples-meta-data/integrated_call_samples_v3.20130502.ALL.panel"
    output:
        "results/data/210305-merged-with-1TGP/pops-to-subset.txt"
    shell:
        """
        cut -f1,2 {input} |\
            grep -E '{oneTGP_pops}' |\
            cut -f1 >{output}
        """


rule file_to_rename_chromosomes:
    output:
        "results/data/210305-merged-with-1TGP/rename-chromosomes.txt"
    shell:
        """
        touch {output}
        for i in {{1..22}}
        do
          echo "chr$i $i" >>{output}
        done
        """


rule oneTGP_get_biallelic_and_subset_pops:
    # subset the population from 1TGP data
    # filter biallelic snps
    # also changes the name of the chromosome from chrN to N
    # to  match the chromosome names in the MXB genomes
    input:
        vcf = path.join(path_to_1TGP, "CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrn}.filtered.shapeit2-duohmm-phased.vcf.gz"),
        pops = "results/data/210305-merged-with-1TGP/pops-to-subset.txt",
        chr_names = "results/data/210305-merged-with-1TGP/rename-chromosomes.txt"
    output:
        temp("results/data/210305-merged-with-1TGP/1TGP-chr{chrn}-biallelic-subset.vcf.gz"),
        temp("results/data/210305-merged-with-1TGP/1TGP-chr{chrn}-biallelic-subset.vcf.gz.tbi")
    # In the bcftools pipeline,
    #I use (at the end) view -c1, this will remove all non variant sites
    shell:
        """
        bcftools view -m2 -M2 -v snps {input.vcf} |\
            bcftools view -S {input.pops} |\
            bcftools annotate --rename-chrs {input.chr_names} |\
            bcftools view -c1 -Oz -o {output[0]}
        bcftools index {output[0]} --tbi
        """


rule merge_50G_with_1TG:
    # merge the 50 genomes with 
    # the 1TGP samples
    # i alse set a variant id for variants
    # with missing id
    input:
        vcf_50g = "results/data/210303-phased/mxb-chr{chrn}-GRCh38-phased.vcf.gz",
        index_50g = "results/data/210303-phased/mxb-chr{chrn}-GRCh38-phased.vcf.gz.tbi",
        vcf_1Tg = "results/data/210305-merged-with-1TGP/1TGP-chr{chrn}-biallelic-subset.vcf.gz",
        index_1Tg = "results/data/210305-merged-with-1TGP/1TGP-chr{chrn}-biallelic-subset.vcf.gz.tbi"
    output:
        temp("results/data/210305-merged-with-1TGP/1TGP_and_50MXB-chr{chrn}-snps.vcf.gz")
    shell:
        """
        bcftools merge --missing-to-ref \
            {input.vcf_1Tg} {input.vcf_50g} |\
            bcftools annotate --set-id \
            +'%CHROM\:%POS\:%REF\:%FIRST_ALT' -Oz -o {output}
        """

# The next rules will add variant annotation information
#with vep
#I use the prefix va_* to mean variant annotation
# To run this rule vep should be installed
# I installed vep in my home directory: ~/ensembl-vep/
# I followed the installation instructions here:
# http://www.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html
# I also downloaded the cache data for GRCh37 and GRCh38

rule va_sort_vcf:
    input:
        "results/data/210305-merged-with-1TGP/1TGP_and_50MXB-chr{chrn}-snps.vcf.gz",
    output:
        temp("results/data/210305-merged-with-1TGP/1TGP_and_50MXB-chr{chrn}-snps-to-annotate.vcf")
    shell:
        """
        bcftools view {input} | bcftools sort >{output}
        """

rule va_vep:
    input:
        "results/data/210305-merged-with-1TGP/1TGP_and_50MXB-chr{chrn}-snps-to-annotate.vcf"
    output:
        "results/data/210305-merged-with-1TGP/1TGP_and_50MXB-chr{chrn}-snps-vep-GRCh38.vcf_summary.html",
        temp("results/data/210305-merged-with-1TGP/1TGP_and_50MXB-chr{chrn}-snps-vep-GRCh38.vcf")
    params:
        assembly="GRCh38",
        species="homo_sapiens"
    log:
        "results/logs/210305-merged-with-1TGP/vep-{chrn}.log"
    shell:
        """
        # NOTE: I am using the path to the vep dir
        ~/ensembl-vep/vep -i {input} \
             --assembly {params.assembly} --cache --vcf \
             --output_file {output[1]} 2>{log}
        """


rule va_compress_and_index:
    input:
        "results/data/210305-merged-with-1TGP/1TGP_and_50MXB-chr{chrn}-snps-vep-GRCh38.vcf"
    output:
        "results/data/210305-merged-with-1TGP/1TGP_and_50MXB-chr{chrn}-snps-vep-GRCh38.vcf.gz",
        "results/data/210305-merged-with-1TGP/1TGP_and_50MXB-chr{chrn}-snps-vep-GRCh38.vcf.gz.tbi"
    shell:
        """
        bcftools view {input} -Oz -o {output[0]}
        bcftools index {output[0]} --tbi
        """

# ****************************************************
# ********** REMOVE MASKED REGIONS ******************* 
# ****************************************************
"""
Why? The masks are useful for population genetic analysis
(such as estimates of mutation rate) that must focus on genomic
regions with very low false positive and false negative rates. 


The StrictMask directory uses a more stringent definition. This definition uses
a narrower band for coverage, requiring that total coverage should be within 50%
of the average, that no more than 0.1% of reads have mapping quality of zero,
and that the average mapping quality for the position should be 56 or greater.

This definition is quite stringent and focuses on the most unique regions of the
genome.
"""


rule masked_bed_regions:
    """
    See the readme file: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/README.accessible_genome_mask.20160622
    This rule gets the bedfile of all passed sites 
    """
    output:
        "results/data/210305-merged-with-1TGP/strict-mask/20160622.allChr.mask.bed" # make this file temporal?, no it is just 16Mb
    shell:
        """
        wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/StrictMask/20160622.allChr.mask.bed
		mv 20160622.allChr.mask.bed {output}
        """

rule masked_chrn:
    """
    Extracts the data for a particular chromosome.
    It also changes the chromosome format name
    from chrXX to XX.
    """
    input:
        "results/data/210305-merged-with-1TGP/strict-mask/20160622.allChr.mask.bed"
    output:
        temp("results/data/210305-merged-with-1TGP/strict-mask/chr{chrn}_mask.bed")
    params:
        chrn="{chrn}"
    shell:
        """
        grep '^chr{params.chrn}' {input} |\
            sed 's/^chr//g' >{output}
        """


rule masked_passed_cites_vcf:
    """
    Filters the masked variants
    from the input vcf
    """
    input:
        mask = "results/data/210305-merged-with-1TGP/strict-mask/chr{chrn}_mask.bed",
        vcf = "results/data/210305-merged-with-1TGP/1TGP_and_50MXB-chr{chrn}-snps-vep-GRCh38.vcf.gz",
        index = "results/data/210305-merged-with-1TGP/1TGP_and_50MXB-chr{chrn}-snps-vep-GRCh38.vcf.gz.tbi"
    output:
        vcf = "results/data/210305-merged-with-1TGP/strict-mask/1TGP_and_50MXB-chr{chrn}-snps-vep-mask-GRCh38.vcf.gz",
        index = "results/data/210305-merged-with-1TGP/strict-mask/1TGP_and_50MXB-chr{chrn}-snps-vep-mask-GRCh38.vcf.gz.tbi"
    shell:
        """
        bcftools view -R {input.mask} {input.vcf} -Oz -o {output[0]}
        bcftools index {output.vcf} --tbi
        """



