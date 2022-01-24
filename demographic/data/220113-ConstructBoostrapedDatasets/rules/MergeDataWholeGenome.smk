'''
Merge bootstrapped data to recover whole genome.
NOTE: I discarded some genome windows because there
were not lof variants there. So this putput is no 100% whole genome,
because it will not include those regions
'''


rule whole_genome_spectrum:
    input:
        expand('data/jSFS/spectrums/spectrum_chunk_{i}_cat_{{varcat}}.pkl.gz', i=BOOSTRAPS)
    output:
        'data/whole-genome/spectrum-cat_{varcat}.pkl.gz'
    params:
        # filename no gz extetions
        pklfile = 'data/whole-genome/spectrum-cat_{varcat}.pkl'
    shell:
        '''
        python scripts/merge-spectrums.py {params.pklfile} {input}
        gzip  {params.pklfile}
        '''

rule whole_genome_mL_noncoding:
    input:
        mls = expand('data/mL-noncoding/mLs/mL_{{ncregion}}_chunk_{i}.txt', i=BOOSTRAPS)
    output:
        'data/whole-genome/mL-{ncregion}.txt'
    shell:
        '''
        cat {input} |\
            cut -d' ' -f2 |\
            awk '{{sum+=$1;}} END{{print sum;}}' >{output}
        '''
        


rule all_whole_genome:
    input:
        expand('data/whole-genome/spectrum-cat_{varcat}.pkl.gz', varcat=CATEGORIES.keys()),
        expand('data/whole-genome/mL-{ncregion}.txt', ncregion=['introns', 'intergenic'])
        
