# Variant annotation mxb-genomes
# To run this rule vep should be installed
# I installed vep in my home directory: ~/ensembl-vep/
# I followed the installation instructions here:
# http://www.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html
# I also downloaded the cache data for GRCh37


rule uncompress_vcf:
    # this rule only generates an uncompressed vcf
    # for vep input, the output is a temporal file
    input:
        "results/data/raw-genomes/mxb-chr{chrn}.vcf.gz"
    output:
        temp("results/data/variant-annotation/mxb-chr{chrn}.vcf")
    shell:
        """
        bcftools view {input} >{output} 
        """


rule variant_annotation:
    input:
        "results/data/variant-annotation/mxb-chr{chrn}.vcf"
    output:
        "results/data/variant-annotation/vep-{chrn}.vcf"
    message: "annotating variants ..."
    log: "results/logs/vep-{chrn}.log"
    params:
        assembly="GRCh37",
        species="homo_sapiens"
    shell:
        """
        # NOTE: I am using the path to the vep dir
        ~/ensembl-vep/vep -i {input} \
            --assembly {params.assembly} --cache --vcf \
             --output_file {output} 2>{log}
        """
  

