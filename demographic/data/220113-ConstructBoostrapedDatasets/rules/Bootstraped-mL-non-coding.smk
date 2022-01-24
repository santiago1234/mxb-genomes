NCREGIONS = ['introns', 'intergenic']

rule nc_get_bootstraped_region:
    # The input file is masked regions
    input:
        regions = '../211128-compute-mL/data/regions/{ncregion}-masked.bed',
        bootstrap = 'data/chunks/chunk_{i}.bed'
    output:
        'data/mL-noncoding/regions/{ncregion}_chunk_{i}.bed'
    params:
        tmpfile = 'tmp-{ncregion}-{i}.bed'
    shell:
        '''
        sed 's/^chr//g' {input.regions} >{params.tmpfile}
        bedtools intersect -a {params.tmpfile} -b {input.bootstrap} >{output}
        rm {params.tmpfile}
        '''


rule nc_get_seq:
    input:
        bed = 'data/mL-noncoding/regions/{ncregion}_chunk_{i}.bed',
        genome = '/data/users/smedina/data-resources/genomes/GRCh38.primary_assembly.genome-nochrprefix.fasta'
    output:
        'data/mL-noncoding/regions/{ncregion}_chunk_{i}.fasta'
    shell:
        '''
		bedtools getfasta -fi {input.genome} -bed {input.bed} >{output}
        '''


rule nc_aggregated_kmers:
    input:
        'data/mL-noncoding/regions/{ncregion}_chunk_{i}.fasta'
    output:
        'data/mL-noncoding/counts/{ncregion}_chunk_{i}.csv'
    threads: 4
    shell:
        '''
        python ../211128-compute-mL/count-3mers-in-seqs.py {input} {output} {threads}
        '''

rule nc_get_mL:
    input:
        'data/mL-noncoding/counts/{ncregion}_chunk_{i}.csv',
        '../211128-compute-mL/data/mutation_rate_methylation_bins.txt'
    output:
        'data/mL-noncoding/mLs/mL_{ncregion}_chunk_{i}.txt'
    shell:
        '''
        python ../211128-compute-mL/mutation-rate.py {input[1]} {input[0]} >{output}
        '''


rule nc_all:
    input:
        expand('data/mL-noncoding/mLs/mL_{ncregion}_chunk_{i}.txt', ncregion=NCREGIONS, i=BOOSTRAPS)
        
        
        

        
        
