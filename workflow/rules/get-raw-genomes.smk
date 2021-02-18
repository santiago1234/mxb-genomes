# copy the genomes and split by chromosome

configfile: "config/config.yaml"

rule split_by_chromosome:
    # to name this vcf files I use the GRCh37 suffix
    # since the built of these genomes is GRCh37
    input:
        config['path_to_raw_genomes']
    output:
        vcf = "results/data/raw-genomes/mxb-chr{chrn}-GRCh37.vcf.gz",
        index = "results/data/raw-genomes/mxb-chr{chrn}-GRCh37.vcf.gz.tbi"
    log: 'results/logs/get-raw-genomes/log-{chrn}-GRCh37.log'
    params:
        chromosome = "{chrn}"
    shell:
        """
        bcftools view -r {params.chromosome} {input} -O b -o {output.vcf}  2> {log}
        bcftools index {output.vcf} --tbi 2>>{log}
        """
