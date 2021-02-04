# Initial quality analysis in mxb genomes

rule all_qc:
    input:
        "results/QC/biallelic-chr{chrn}.vcf.gz",
        "results/QC/tmp-dir/chr{chrn}-seqdepth.csv",
        "results/QC/tmp-dir/chr{chrn}-vars-per-sample.txt"

rule get_biallelic_snps:
    # for QC we only analyze biallelic loci
    input:
        "results/data/raw-genomes/mxb-chr{chrn}.vcf.gz"
    output:
        "results/QC/biallelic-chr{chrn}.vcf"
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
        "results/QC/tmp-dir/chr{chrn}-vars-per-sample.txt"
    message: "Counting variants in {input}"
    shell:
        """
        mxb_samples=`bcftools query -l {input}`
        touch {output}
        for sample in $mxb_samples
        do
            n_vars=`bcftools view -c1 -H -s $sample {input} |cut -f1 |uniq -c`
            echo $sample $n_vars >>{output}
        done
        """

rule aggregate_qc_data:
    # aggregate the qc data for each chromosome
    input:
        seq_deps = expand("results/QC/tmp-dir/chr{chrn}-seqdepth.csv", chrn=[20, 21, 22]),
        n_vars = expand("results/QC/tmp-dir/chr{chrn}-vars-per-sample.txt", chrn=[20, 21, 22])
      
