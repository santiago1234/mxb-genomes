# Local Ancestry with XGMix
#
# For reference populations:
# AFR: YRI (1TGP)
# EUR: IBR + GBR (1TGP)
# EAS: CHB (1TGP)
# NAT: MXB
#
# Query pop
# MXL + PEL + CLM + PUR: AMR (1TGP)

oneTGP_pops = config['oneTGP_pops']


rule download_XGMix:
    # this rules dowloads the software (git repo)
    # to run XGMix
    output:
        directory("workflow/scripts/XGMix-master")
    shell:
        """
        wget https://github.com/AI-sandbox/XGMix/archive/master.zip
        unzip master.zip
        rm master.zip
        mv XGMix-master workflow/scripts/
        """


rule sample_map_and_query_pops:
    # sample_map_file is a sample map file matching reference samples to their respective reference populations
    # query pops is a file with the query populations
    # the script 01-populations.R will degine the reference and query haplotypes
    input:
        oneT_meta = "resources/1TGP-samples-meta-data/igsr-1000genomes.tsv",
        mxb_meta = "resources/genomes-metadata/50Genomes_info.txt"
    output:
        sample_map_file = "results/data/210308-local-ancestry/input-data/sample_map_file.txt",
        query_pops = "results/data/210308-local-ancestry/input-data/query-pops.txt"
    params:
        oneTGP_pops = oneTGP_pops
    conda: "../envs/renv-yaml"
    script:
        "../scripts/local-ancestry/01-populations.R"


