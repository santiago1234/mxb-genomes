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
        
