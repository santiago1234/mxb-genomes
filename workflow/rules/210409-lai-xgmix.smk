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
    script:
        '../scripts/local-ancestry/define_query_and_reference_populations.py'
