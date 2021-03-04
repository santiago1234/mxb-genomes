# PHASING 50 MXB-genomes
# note: phasing is done in the GRCh38 genome build


rule download_genetic_maps:
    # download genetic maps from github repo:
    # https://github.com/joepickrell/1000-genomes-genetic-maps
    output:
        expand("resources/genetic-maps/chr{chrn}.b38.gmap", chrn=CHROMS)
    shell:
        """
        wget https://github.com/odelaneau/shapeit4/raw/master/maps/genetic_maps.b38.tar.gz
        tar -xzf genetic_maps.b38.tar.gz
        gunzip chr*.gz
        mkdir -p resources/genetic-maps/
        mv *gmap resources/genetic-maps/
        rm genetic_maps.b38.tar.gz
        """


rule biallelic_snps_to_phase:
    # SHAPEIT2 does not handle multiallelic variant phasing
    # this rule retrieves biallelic snp
    # and removes duplicated positions. For example chromosome 22 has
    # ~22 duplicated positions. This is a minimal number but shapeit is not
    # able to run with duplicated positions
    input:
        "results/data/raw-genomes/mxb-chr{chrn}-GRCh38.vcf.gz"
    output:
        temp("results/data/210303-phased/mxb-unphased-chr{chrn}-GRCh38.vcf.gz")
    shell:
        """
        bcftools view -m2 -M2 -v snps {input} |\
            bcftools norm -d none -Oz -o {output}
        """


rule phase:
    input:
        vcf = "results/data/210303-phased/mxb-unphased-chr{chrn}-GRCh38.vcf.gz",
        genetic_map = "resources/genetic-maps/chr{chrn}.b38.gmap"
    output:
        temp("results/data/210303-phased/mxb-chr{chrn}-GRCh38-phased.sample"),
        temp("results/data/210303-phased/mxb-chr{chrn}-GRCh38-phased.haps")
    log:
        "results/logs/phasing/mxb-chrn{chrn}"
    params:
        out_basename = "results/data/210303-phased/mxb-chr{chrn}-GRCh38-phased"
    threads:
        5
    shell:
        """
        shapeit --input-vcf {input.vcf} -O {params.out_basename} \
            -T {threads} -M {input.genetic_map} --output-log {log}
        """


rule convert_haps_to_vcf:
    input:
        "results/data/210303-phased/mxb-chr{chrn}-GRCh38-phased.sample",
        "results/data/210303-phased/mxb-chr{chrn}-GRCh38-phased.haps"
    output:
        "results/data/210303-phased/mxb-chr{chrn}-GRCh38-phased.vcf.gz"
    params:
        basename = "results/data/210303-phased/mxb-chr{chrn}-GRCh38-phased",
        out_vcf = "results/data/210303-phased/mxb-chr{chrn}-GRCh38-phased.vcf"
    shell:
        """
        shapeit -convert --input-haps {params.basename} \
            --output-vcf {params.out_vcf}
        bcftools view {params.out_vcf} -Oz -o {output}
        rm {params.out_vcf}
        """
