rule cr_exons_in_region:
    input:
        # NOTE: We are passing the non-masked coding regions of the genome
        exons = '../220423-Coding-mL-Whole-Genome-test/data/coding-regions-sorted.bed',
        regions = 'data/samples/region_region_{i}.bed'
    output:
        'data/samples/region_exons_{i}.bed'
    shell:
        '''
        bedtools intersect -a {input.exons} -b {input.regions} |\
            sort -k 1,1n -k2,2n |\
            bedtools merge -i -  >{output}
        '''

rule cr_all_possible_variants_in_reg:
    input:
        bed = 'data/samples/region_exons_{i}.bed',
        genome = '/data/users/smedina/data-resources/genomes/GRCh38.primary_assembly.genome-nochrprefix.fasta'
    output:
        'data/ml-coding/variants/region_{i}.bed.gz'
    params:
         outbed = 'tmp-vars-chunk_{i}.bed'
    shell:
        '''
        python ../211231-mL-coding/scripts/all-snps-in-regions.py {input.bed} {input.genome} {params.outbed}
        sort {params.outbed} |\
            uniq |\
            sort -k 2,2n |\
            gzip >{output}
        rm -f {params.outbed}
        '''


rule cr_predict_variants:
    input:
        'data/ml-coding/variants/region_{i}.bed.gz'
    output:
        'data/ml-coding/vep/region_{i}.txt.gz',
        'data/ml-coding/vep/region_{i}.txt.gz_warnings.txt',
        'data/ml-coding/vep/region_{i}.txt.gz_summary.html'
    threads: 5
    shell:
        '''
        ~/ensembl-vep/vep -i {input} --cache  \
             --fork {threads} --verbose \
            --assembly GRCh38 --tab --output_file {output[0]} \
            --compress_output gzip --coding_only  --fields \
            'Uploaded_variation,Location,Allele,Gene,Feature_type,Consequence,Codons'
        '''


rule cr_process_vep_output:
    input:
        'data/ml-coding/vep/region_{i}.txt.gz',
    output:
        'data/ml-coding/vep/region-UNIQ_{i}.txt.gz',
    shell:
        '''
        python ../211231-mL-coding/scripts/process-vep-output.py {input} {output}
        '''

rule cr_compute_ml:
    input:
        vep = 'data/ml-coding/vep/region-UNIQ_{i}.txt.gz',
        mus = '../211128-compute-mL/data/mutation_rate_methylation_bins.txt' # mutation rates
    output:
        'data/samples/region_mlcoding_{i}.csv'
    shell:
        '''
        python ../211231-mL-coding/scripts/compute-mL-from-VEP.py {input.vep} {input.mus} {wildcards.i} {output}
        '''

        
