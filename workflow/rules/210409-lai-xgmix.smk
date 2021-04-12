# LAI with XGMix

# Here, I run XGMix indpendently twice
# 3pops:
#   A three population model, using the following reference panels
#       - EUR: IBR + GBR
#       - NAT: 50 MXB + 1TGP (PEL and MXL that were shown to be representative for NAT ancestry)
#       - AFR: YRI
# 4pops:
#   A four population model, using the following reference panels
#       - EUR: IBR + GBR
#       - NAT: 50 MXB + 1TGP (PEL and MXL that were shown to be representative for NAT ancestry)
#       - AFR: YRI
#       - EAS: CHB


rule lai_ref_and_query_pops:
    input:
     # reference_map_file is a sample map file matching reference samples to their respective reference populations
     # query_map_file is a file with the query populations
        oneT_meta = "resources/1TGP-samples-meta-data/integrated_call_samples_v3.20130502.ALL.panel",
        mxb_meta = "resources/genomes-metadata/50Genomes_info.txt",
        oneT_native_american = "resources/1TGP-samples-meta-data/native-american.txt"
    output:
        sample_map_file = "results/data/210409-local-ancestry/{npops}-pops/input-data/reference_map_file.txt",
        query_pops = "results/data/210409-local-ancestry/{npops}-pops/input-data/query_map_file.txt"
    params:
        npops = '{npops}'
    conda:
        "../envs/scipy.yaml"
    script:
        '../scripts/local-ancestry/define_query_and_reference_populations.py'


rule lai_pop_list_to_subset:
# generates a list of populations
# to subset from a vcf fiel
    input:
        "results/data/210409-local-ancestry/{npops}-pops/input-data/{haplo}_map_file.txt"
    output:
        "results/data/210409-local-ancestry/{npops}-pops/input-data/{haplo}_pop-list.txt"
    shell:
        """
        cut -f1 {input} |\
            grep -v "Sample" >{output} #to remove the header
        """


rule lai_subset_haplotypes:
    input:
        vcf = "results/data/210305-merged-with-1TGP/1TGP_and_50MXB-chr{chrn}-snps-vep-GRCh38.vcf.gz",
        pop_list = "results/data/210409-local-ancestry/{npops}-pops/input-data/{haplo}_pop-list.txt"
    output:
        temp("results/data/210409-local-ancestry/{npops}-pops/input-data/{haplo}-chr{chrn}.vcf.gz")
    shell:
        """
        bcftools view -S {input.pop_list} {input.vcf} |\
            bcftools view -c1 -Oz -o {output}
        """

rule lai_genetic_map:
#Xgmix expects a tsv file that has the following columns:
#chm pos pos_cm
#here I just reformat the genetic map that is available
#to be used with XGMix
    input:
        "resources/genetic-maps/chr{chrn}.b38.gmap"
    output:
        temp("results/data/210409-local-ancestry/chr{chrn}.b38.gmap.txt")
    shell:
        """
        awk -v OFS='\t' '{{print $2, $1, $3}}' {input} |\
          grep -v '^chr' >{output} #XGMix does not expect a header
        """


rule lai_xgmix_train:
    input:
        sample_map_file = "results/data/210409-local-ancestry/{npops}-pops/input-data/reference_map_file.txt",
        query_file = "results/data/210409-local-ancestry/{npops}-pops/input-data/query-chr{chrn}.vcf.gz",
        reference_file = "results/data/210409-local-ancestry/{npops}-pops/input-data/reference-chr{chrn}.vcf.gz",
        genetic_map = "results/data/210409-local-ancestry/chr{chrn}.b38.gmap.txt"
    params:
        phase="True",
        chr_nr="{chrn}",
        output_basename="results/data/210409-local-ancestry/{npops}-pops/mdl-{chrn}"
    conda: "../envs/xgmix.yaml"
    log: "results/logs/210409-xgmix/{npops}-{chrn}-xgmix.log"
    threads: 20
    output:
        temp(directory("results/data/210409-local-ancestry/{npops}-pops/mdl-{chrn}"))
    shell:
        """
        echo "RUNING XGmix !!!!!"
        python workflow/scripts/XGMix-master/XGMIX.py {input.query_file} \
            {input.genetic_map} {params.output_basename} \
            {params.chr_nr} {params.phase} {input.reference_file} \
            {input.sample_map_file} 2>{log}
        """


rule lai_collect_xgmix_output:
# this rule simply moves the important output generated by XMGix
# to another directory. This is to save disk space, since each model
# gnerated ~100G of data!. The XGMix output dir will be deleted.
    input:
        "results/data/210409-local-ancestry/{npops}-pops/mdl-{chrn}"
    output:
        mdl = "results/data/210409-local-ancestry/{npops}-pops/trained-models/model_chm_{chrn}.pkl",
        pred_msp = "results/data/210409-local-ancestry/{npops}-pops/predictions/mdl-{chrn}.msp.tsv",
        pred_fb = "results/data/210409-local-ancestry/{npops}-pops/predictions/mdl-{chrn}.fb.tsv",
        rephased_vcf = "results/data/210409-local-ancestry/{npops}-pops/rephased-query-with-lai/query-chr{chrn}.vcf.gz"
    params:
        #these are mainly output files generated by xgmixi
        xmdl = "results/data/210409-local-ancestry/{npops}-pops/mdl-{chrn}/models/model_chm_{chrn}/model_chm_{chrn}.pkl",
        xpred_msp = "results/data/210409-local-ancestry/{npops}-pops/mdl-{chrn}/mdl-{chrn}.msp.tsv",
        xpred_fb = "results/data/210409-local-ancestry/{npops}-pops/mdl-{chrn}/mdl-{chrn}.fb.tsv",
        xrephase_vcf = "results/data/210409-local-ancestry/{npops}-pops/mdl-{chrn}/query_file_phased.vcf"
    shell:
        """
        mv {params.xmdl} {output.mdl}
        mv {params.xpred_msp} {output.pred_msp}
        mv {params.xpred_fb} {output.pred_fb}
        bcftools view {params.xrephase_vcf} -Oz -o {output.rephased_vcf}
        """