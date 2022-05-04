rule sr_make_genome_windows:
    input:
        '../../data/220113-ConstructBoostrapedDatasets/data/human-autosomes.genome'
    output:
        'data/chunks/genome-in-1MB-non-overlaping-windows.bed'
    shell:
        '''
        # bedtools complains for having column names
        tail -n +2 {input} |\
            bedtools makewindows -g - -w 100000 >{output}
        '''


rule sr_region:
    input:
        "data/chunks/genome-in-1MB-non-overlaping-windows.bed"
    output:
        temp(expand('data/chunks/chunk_{i}', i=N_CHUNKS))
    shell:
        """
        python scripts/define-windows.py {input}
        """


rule sr_merge_region:
    # Merge, this is just for a nicer format in output files.
    input:
        'data/chunks/chunk_{i}'
    output:
        'data/regions/region_{i}.bed'
    shell:
        '''
        bedtools merge -i {input} >{output}
        '''


rule sr_all:
    input:
        expand('data/regions/region_{i}.bed', i=N_CHUNKS)
