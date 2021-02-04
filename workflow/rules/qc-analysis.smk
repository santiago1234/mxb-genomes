# Initial quality analysis in mxb genomes

rule all_qc:
    input:
        "results/QC/biallelic-chr{chrn}.vcf.gz"

rule get_biallelic_snps:
    # for QC we only analyze biallelic loci
    input:
        "results/data/raw-genomes/mxb-chr{chrn}.vcf.gz"
    output:
        temp("results/QC/biallelic-chr{chrn}.vcf")
    shell:
        """
        bcftools view -m2 -M2 -v snps {input} > {output}
        """


rule sequence_depth:
    # generates tables with sequence depth, etc.
    # for qc plots
    input:
        "results/QC/biallelic-chr{chrn}.vcf"
    output:
        "results/QC/tmp-dir/chr{chrn}-seqdepth.csv"
    conda:
        "../envs/renv.yaml"
    params:
        # the fraction of variants to get from each chromosome
        fraction=0.1
    script:
        "../scripts/qc/get-sequence-depth.R"


rule count_variants_per_sample:
    # count the variants per sample
    # see: https://www.biostars.org/p/336206/
    input:
        "results/QC/biallelic-chr{chrn}.vcf",
    output:
        temp("results/QC/tmp-dir/chr{chrn}-vars-per-sample.txt")
    message: "Counting variants in {input}"
    shell:
        """
        # list the samples from the vcf file
        # the output format is indvidual-id n_variants chr
        mxb_samples=`bcftools query -l {input}`
        # create empty output file
        touch {output}
        for sample in $mxb_samples
        do
            n_vars=`bcftools view -c1 -H -s $sample {input} |cut -f1 |uniq -c`
            echo $sample $n_vars >>{output}
        done
        """

