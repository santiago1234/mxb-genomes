rule mnc_noncoding_region:
    input:
        # These regions are not masked
        introns = '../211128-compute-mL/data/regions/introns.bed',
        intergenic = '../211128-compute-mL/data/regions/intergenic.bed'
    output:
        'data/ml-noncoding/non-coding-regions.bed'
    shell:
        '''
        cat {input} |\
            sed 's/^chr//g' >non-coding-tmp.bed

        sort -k 1,1n -k2,2n non-coding-tmp.bed >non-coding-tmp-sorted.bed

        bedtools merge -i non-coding-tmp-sorted.bed >{output}

        rm non-coding-tmp-sorted.bed non-coding-tmp.bed
        '''


rule mnc_get_reg:
    input:
        nc = 'data/ml-noncoding/non-coding-regions.bed',
        region = 'data/samples/region_region_{i}.bed'
    output:
        'data/samples/region_intronANDinterg_{i}.bed'
    shell:
        '''
        bedtools intersect -a {input.nc} -b {input.region} >{output}
        '''

rule mnc_get_seq:
    input:
        bed = 'data/samples/region_intronANDinterg_{i}.bed',
        genome = '/data/users/smedina/data-resources/genomes/GRCh38.primary_assembly.genome-nochrprefix.fasta'
    output:
        'data/ml-noncoding/region_{i}.fasta'
    shell:
        '''
		bedtools getfasta -fi {input.genome} -bed {input.bed} >{output}
        '''

rule mnc_aggregated_kmers:
    input:
        'data/ml-noncoding/region_{i}.fasta'
    output:
        'data/ml-noncoding/counts/region_{i}.csv'
    threads: 4
    shell:
        '''
        python ../211128-compute-mL/count-3mers-in-seqs.py {input} {output} {threads}
        '''

rule mnc_get_mL:
    input:
        'data/ml-noncoding/counts/region_{i}.csv',
        '../211128-compute-mL/data/mutation_rate_methylation_bins.txt'
    output:
        'data/samples/region_mlnoncoding_{i}.txt'
    shell:
        '''
        python ../211128-compute-mL/mutation-rate.py {input[1]} {input[0]} >{output}
        '''
