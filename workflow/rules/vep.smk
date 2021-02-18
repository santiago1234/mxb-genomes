# Variant annotation mxb-genomes
# To run this rule vep should be installed
# I installed vep in my home directory: ~/ensembl-vep/
# I followed the installation instructions here:
# http://www.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html
# I also downloaded the cache data for GRCh37 and GRCh38

rule sort_vcf:
    # vep needs the vcf file to be sorted
    input:
        "results/data/raw-genomes/mxb-chr{chrn}-GRCh{build}.vcf.gz"
    output:
        temp("results/data/variant-annotation/sorted-chr{chrn}-GRCh{build}.vcf")
    shell:
        """
        bcftools view {input} | bcftools sort >{output}
        """


rule variant_annotation:
    input:
        "results/data/variant-annotation/sorted-chr{chrn}-GRCh{build}.vcf"
    output:
        "results/data/variant-annotation/mxb-chr{chrn}-GRCh{build}.vcf_summary.html",
        temp("results/data/variant-annotation/mxb-chr{chrn}-GRCh{build}.vcf")
    message: "annotating variants ..."
    log: "results/logs/vep/mxb-chr{chrn}-GRCh{build}.log"
    params:
        assembly="GRCh{build}",
        species="homo_sapiens",
    shell:
        """
        # NOTE: I am using the path to the vep dir
        ~/ensembl-vep/vep -i {input} \
            --assembly {params.assembly} --cache --vcf \
             --output_file {output[1]} 2>{log}
        """


rule compress_annotated_vcf:
    input:
        "results/data/variant-annotation/mxb-chr{chrn}-GRCh{build}.vcf"
    output:
        "results/data/variant-annotation/mxb-chr{chrn}-GRCh{build}.vcf.gz"
    shell:
        """
        bcftools view {input} -O b -o {output}
        """
  

