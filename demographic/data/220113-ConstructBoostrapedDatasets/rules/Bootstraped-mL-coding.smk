CATEGORIES = ['LOF', 'missense', 'synonymous']

rule cod_get_boostraped_region:
   """
   NOTE: The inpute is the file with the exons after
   applying the MASK. See the Snakemake file in: ../211231-mL-coding/
   
   Then I intersect this file with the bootstrap regions
   """
   input:
       exons_masked = '../211231-mL-coding/data/regions/exones-mask.bed',
       bootstrap = 'data/chunks/chunk_{i}.bed'
   output:
       'data/mL-coding/regions/exons_chunk_{i}.bed'
   params:
       # I need to add the chr prefix to this file
       # so the chromosome names are consistent
       tmpfile = 'tmp-chunk-chrprefix-{i}.bed'
   shell:
       """
       # First add the chr prefix and make tmp file
       cat {input.bootstrap} |awk '{{print "chr" $0}}' >{params.tmpfile}
       bedtools intersect -a {input.exons_masked} -b {params.tmpfile} >{output}
       rm {params.tmpfile}
       """


rule cod_get_all_possible_variants_in_Regions:
    """
    This rule generates every possible variant so we can run VEP.
    """
    input:
       bed = 'data/mL-coding/regions/exons_chunk_{i}.bed',
       genome = '/data/users/smedina/data-resources/genomes/GRCh38.primary_assembly.genome-nochrprefix.fasta'
    output:
        'data/mL-coding/tmp/variants/chunk_{i}.bed.gz'
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

rule cod_predict_variants:
    input:
        'data/mL-coding/tmp/variants/chunk_{i}.bed.gz'
    output:
        'data/mL-coding/tmp/vep/chunk_{i}.txt.gz',
        'data/mL-coding/tmp/vep/chunk_{i}.txt.gz_warnings.txt',
        'data/mL-coding/tmp/vep/chunk_{i}.txt.gz_summary.html'
    threads: 5
    shell:
        '''
        ~/ensembl-vep/vep -i {input} --cache  \
             --fork {threads} --verbose \
            --assembly GRCh38 --tab --output_file {output[0]} \
            --compress_output gzip --coding_only  --fields \
            'Uploaded_variation,Location,Allele,Gene,Feature_type,Consequence,Codons'
        '''
        

rule cod_process_vep_output:
    input:
        'data/mL-coding/tmp/vep/chunk_{i}.txt.gz',
    output:
        'data/mL-coding/tmp/vep/chunk-UNIQ_{i}.txt.gz',
    shell:
        '''
        python ../211231-mL-coding/scripts/process-vep-output.py {input} {output}
        '''


rule cod_compute_mL:
    input:
        vep = 'data/mL-coding/tmp/vep/chunk-UNIQ_{i}.txt.gz',
        mus = '../211128-compute-mL/data/mutation_rate_methylation_bins.txt' # mutation rates
    output:
        'data/mL-coding/tmp/mLs/mLs-csqs_{i}.csv'
    shell:
        '''
        python ../211231-mL-coding/scripts/compute-mL-from-VEP.py {input.vep} {input.mus} {wildcards.i} {output}
        '''

rule cod_synonymous_mL:
    """
    This rule simply extract the synonymous mL from the table above.
    """
    input:
        'data/mL-coding/tmp/mLs/mLs-csqs_{i}.csv'
    output:
        'data/mL-coding/mLs/mL_{cat}_chunk_{i}.txt'
    shell:
        '''
        grep '{wildcards.cat}' {input} |\
            tr ',' ' ' |\
            cut -d' ' -f1,2 |\
            sed  's/{wildcards.cat}/mL:/' >{output}
        '''

rule cod_all:
    input:
        expand('data/mL-coding/mLs/mL_{cat}_chunk_{i}.txt', cat=CATEGORIES, i=BOOSTRAPS)
