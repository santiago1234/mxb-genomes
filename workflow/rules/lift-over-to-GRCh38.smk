# Lift over vcf files to GRCh38
# I use the programs crossmap
# documentation: http://crossmap.sourceforge.net/#convert-vcf-format-files
# I create a conda enviroment: envs/crossmap.yaml 


rule download_GRCh38_genome:
    # https://www.biostars.org/p/271395/
    # I donwload the genome from genecode
    # https://www.gencodegenes.org/human/
    output:
        temp("results/data/lifted_to_GRCh38/GRCh38.primary_assembly.genome.fa")
    shell:
        """
        wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz
        gunzip GRCh38.primary_assembly.genome.fa.gz
        mv GRCh38.primary_assembly.genome.fa results/data/lifted_to_GRCh38/
        """


rule download_ucsc_chain_hg19ToHg38:
    # http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
    output:
        temp("results/data/lifted_to_GRCh38/hg19ToHg38.over.chain.gz")
    shell:
        """
        wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
        mv hg19ToHg38.over.chain.gz results/data/lifted_to_GRCh38/
        """


rule lift_over_to_GRCh38:
    input:
        vcf = config['path_to_raw_genomes'],
        chain_file = "results/data/lifted_to_GRCh38/hg19ToHg38.over.chain.gz",
        GRCh38_genome = "results/data/lifted_to_GRCh38/GRCh38.primary_assembly.genome.fa"
    output:
        "results/data/lifted_to_GRCh38/mxb-ALL-GRCh38.vcf.gz",
        "results/data/lifted_to_GRCh38/mxb-ALL-GRCh38.vcf.unmap"
    log: "results/logs/lift-over-to-GRCh38/log.log"
    conda: "../envs/crossmap.yaml"
    params:
        # a temporal file to uncompress the vcf
        vcf_uncompressed = "results/data/lifted_to_GRCh38/mxb-ALL-GRCh37.vcf",
        out_vcf = "results/data/lifted_to_GRCh38/mxb-ALL-GRCh38.vcf"
    shell:
        """
        bcftools view {input.vcf} >{params.vcf_uncompressed}
        CrossMap.py vcf {input.chain_file} {params.vcf_uncompressed} \
            {input.GRCh38_genome} {params.out_vcf} 2>{log}
        rm -f {params.vcf_uncompressed}
        bcftools view {params.out_vcf} | bcftools sort -O b -o {output[0]}
        rm -f {params.out_vcf}
        """


rule index_lifted_vcf:
    input:
        "results/data/lifted_to_GRCh38/mxb-ALL-GRCh38.vcf.gz"
    output: "results/data/lifted_to_GRCh38/mxb-ALL-GRCh38.vcf.gz.tbi"
    shell:
        """
        bcftools index -t {input}
        """


rule lifted_split_by_chromosome:
    # this rule merges the unmapped variant intp on vcf
    input:
        "results/data/lifted_to_GRCh38/mxb-ALL-GRCh38.vcf.gz",
        "results/data/lifted_to_GRCh38/mxb-ALL-GRCh38.vcf.gz.tbi"
    output:
        "results/data/raw-genomes/mxb-chr{chrn}-GRCh38.vcf.gz"
    params:
            chromosome = "{chrn}"
    shell:
        """
        bcftools view -r {params.chromosome} {input[0]} | bcftools sort -O b -o {output}
        """
