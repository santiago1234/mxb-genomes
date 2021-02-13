#
# Initial quality analysis in mxb genomes

# workflow parameters

FRACTION = 0.1

rule all_qc:
    input:
        "results/plots/qc/vars_per_genome.png",
        "results/plots/qc/depth_per_sample.png",
        "results/plots/qc/depth_in_chr22.png",
        "results/plots/qc/missing_data_by_ind.png",
        "results/plots/qc/missing_data_by_var.png",
        "results/plots/qc/pca_comp1_2_3.png",
        "results/plots/qc/pca_withEthnicity_Geography.png"
    output:
        "results/QC/note.txt"
    shell:
        """
        # this file is an indicator for pipeline completion
        date_today=$(date)
        echo "QC completed on $(date_today)" >{output}
        """

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
        temp("results/QC/tmp-dir/chr{chrn}-seqdepth.csv")
    conda:
        "../envs/renv.yaml"
    params:
        # the fraction of variants to get from each chromosome
        fraction=FRACTION
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
        seq_deps = expand("results/QC/tmp-dir/chr{chrn}-seqdepth.csv", chrn=CHROMS),
        n_vars = expand("results/QC/tmp-dir/chr{chrn}-vars-per-sample.txt", chrn=CHROMS)
    output:
        vars_per_genome = "results/QC/nvars_per_genome.csv",
        seqs_deps = "results/QC/sequence-depth.csv"
    conda:
        "../envs/renv.yaml"
    script:
        "../scripts/qc/aggregate_qc_data.R"


rule qc_plots:
    input:
        vars_per_genome = "results/QC/nvars_per_genome.csv",
        seqs_deps = "results/QC/sequence-depth.csv",
    output:
        vars_per_genome_plt = "results/plots/qc/vars_per_genome.png",
        depth_per_sample_plt = "results/plots/qc/depth_per_sample.png",
        depth_in_chr22_plt = "results/plots/qc/depth_in_chr22.png",
        miss_ind_plt = "results/plots/qc/missing_data_by_ind.png",
        miss_var_plt = "results/plots/qc/missing_data_by_var.png"
    message: "generating plots QC"
    log: "results/logs/qc-plot.log"
    conda:
        "../envs/renv.yaml"
    script:
        "../scripts/qc/plots.R"


# PCA plot ----------------------------------------------------------------
# the following rules are used to generate a pca plot


rule convert2plink_and_filter:
    input:
        expand("results/QC/biallelic-chr{chrn}.vcf", chrn=CHROMS)
    output:
        temp(multiext("results/QC/biallelic-ALL", ".fam", ".bed", ".bim"))
    params:
        temp_vcf_all_chrn = "results/QC/biallelic-ALL.vcf",
        plink_out_file = "results/QC/biallelic-ALL",
        # filters
        maf = 0.05,
        hwe = 0.001
    log: "results/logs/QC/convert2plink.log"
    shell:
        """
        # I use bcftools to concat the chromosomes into a
        # single vcf, then I use the annotate command
        # to and id to the variant, see:
        # https://github.com/samtools/bcftools/issues/178
        # the reason that id add a variant id is that
        # it is usefull for plink
        bcftools concat {input} |\
            bcftools annotate \
            --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
            -o {params.temp_vcf_all_chrn} 2>{log}

        plink --vcf {params.temp_vcf_all_chrn} \
            --keep-allele-order \
            --maf {params.maf} \
            --hwe {params.hwe} \
            --make-bed --out {params.plink_out_file} 2>>{log}
        rm {params.temp_vcf_all_chrn}  results/QC/biallelic-ALL.log results/QC/biallelic-ALL.nosex
        """


rule snps_in_linkage_eq:
    input:
        multiext("results/QC/biallelic-ALL", ".fam", ".bed", ".bim")
    output:
        temp(multiext("results/QC/ldpruned_snplist", ".prune.in", ".prune.out")),
        temp(multiext("results/QC/ldpruned_data", ".fam", ".bed", ".bim"))
    shell:
        """
        plink --bfile results/QC/biallelic-ALL \
            --indep-pairwise 200kb 1 0.5 \
            --out results/QC/ldpruned_snplist

        plink --bfile results/QC/biallelic-ALL \
            --extract results/QC/ldpruned_snplist.prune.in \
            --make-bed \
            --out results/QC/ldpruned_data
        """


rule qc_pca:
    input:
        multiext("results/QC/ldpruned_data", ".fam", ".bed", ".bim")
    output:
        multiext("results/QC/pca_results", ".eigenvec", ".eigenval")
    shell:
        """
        plink --bfile results/QC/ldpruned_data  --pca 20 --out results/QC/pca_results
        rm -f results/QC/*.log results/QC/*.nosex
        """


rule qc_pca_plot:
    input:
        multiext("results/QC/pca_results", ".eigenvec", ".eigenval"),
        "resources/genomes-metadata/50Genomes_info.txt"
    output:
        "results/plots/qc/pca_comp1_2_3.png",
        "results/plots/qc/pca_withEthnicity_Geography.png"
    conda:
        "../envs/renv.yaml"
    script:
        "../scripts/qc/plot_pca.R"


