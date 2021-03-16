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
    # parameters for runing XGMix can be configured in the
    # config file inside XGMix-master
    # I only manually change the follwing parameters:
    #  rm_simulated_data = True
    output:
        directory("workflow/scripts/XGMix-master")
    shell:
        """
        wget https://github.com/AI-sandbox/XGMix/archive/master.zip
        unzip master.zip
        rm master.zip
        mv XGMix-master workflow/scripts/
        """


rule reference_and_query_pops:
    # reference_map_file is a sample map file matching reference samples to their respective reference populations
    # query_map_file is a file with the query populations
    # the script 01-populations.R will degine the reference and query haplotypes
    input:
        oneT_meta = "resources/1TGP-samples-meta-data/igsr-1000genomes.tsv",
        mxb_meta = "resources/genomes-metadata/50Genomes_info.txt"
    output:
        sample_map_file = "results/data/210308-local-ancestry/input-data/reference_map_file.txt",
        query_pops = "results/data/210308-local-ancestry/input-data/query_map_file.txt"
    params:
        oneTGP_pops = oneTGP_pops
    conda: "../envs/renv.yaml"
    script:
        "../scripts/local-ancestry/01-populations.R"


rule pop_list_to_subset:
    input:
        "results/data/210308-local-ancestry/input-data/{haplo}_map_file.txt"
    output:
        temp("results/data/210308-local-ancestry/input-data/{haplo}_pop-list.txt")
    shell:
        """
        cut -f1 {input} |\
            grep -v "Sample" >{output} #to remove the header
        """


rule subset_haplotypes:
    input:
        vcf = "results/data/210305-merged-with-1TGP/1TGP_and_50MXB-chr{chrn}-snps-GRCh38.vcf.gz",
        pop_list = "results/data/210308-local-ancestry/input-data/{haplo}_pop-list.txt"
    output:
        temp("results/data/210308-local-ancestry/input-data/{haplo}-chr{chrn}.vcf.gz")
    params:
        # temporal file to save ref pop list
        pop_list_tmp = "results/data/210308-local-ancestry/sample_pop_list_tmp"
    shell:
        """
        bcftools view -S {input.pop_list} {input.vcf} -Oz -o {output}
        """


rule xgmix_genetic_map:
    #Xgmix expects a tsv file that has the following columns:
    #chm pos pos_cm
    #here I just reformat the genetic map that is availanle
    #to be used with XGMix
    input:
        "resources/genetic-maps/chr{chrn}.b38.gmap"
    output:
        temp("results/data/210308-local-ancestry/input-data/chr{chrn}.b38.gmap.txt")
    shell:
        """
        awk -v OFS='\t' '{{print $2, $1, $3}}' {input} |\
            grep -v '^chr' >{output} #XGMix does not expect a header
        """


rule xgmix_train:
    input:
        sample_map_file = "results/data/210308-local-ancestry/input-data/reference_map_file.txt",
        query_file = "results/data/210308-local-ancestry/input-data/query-chr{chrn}.vcf.gz",
        reference_file = "results/data/210308-local-ancestry/input-data/reference-chr{chrn}.vcf.gz",
        genetic_map = "results/data/210308-local-ancestry/input-data/chr{chrn}.b38.gmap.txt"
    params:
        phase="False",
        chr_nr="{chrn}",
        output_basename="results/data/210308-local-ancestry/mdl-{chrn}"
    conda: "../envs/xgmix.yaml"
    log: "results/logs/210308-xgmix/{chrn}-xgmix.log"
    threads: 20
    output:
        directory("results/data/210308-local-ancestry/mdl-{chrn}")
    shell:
        """
        echo "RUNING XGmix !!!!!"
        python workflow/scripts/XGMix-master/XGMIX.py {input.query_file} \
            {input.genetic_map} {params.output_basename} \
            {params.chr_nr} {params.phase} {input.reference_file} \
            {input.sample_map_file} 2>{log}
        """

# rule collect_xgmix_output:


