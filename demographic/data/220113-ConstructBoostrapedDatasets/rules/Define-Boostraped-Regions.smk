DEFINED_BOOSTRAPS = list(range(1, N_BOOSTRAPS + 1))

rule make_genome_windows:
    """
    The human autosome genome length is 2875001522
    """
    input:
        "data/human-autosomes.genome"
    output:
        "data/chunks/genome-in-1MB-non-overlaping-windows.bed"
    shell:
        """
        bedtools makewindows -g {input} -w 100000 >{output}
        """


rule bootstrap_windows:
    input:
        "data/chunks/genome-in-1MB-non-overlaping-windows.bed"
    output:
        temp(expand('data/chunks/chunk_{i}', i=DEFINED_BOOSTRAPS))
    shell:
        """
        python scripts/define-windows.py {input}
        """


rule merge:
    # Merge, this is just for a nicer format in output files.
    input:
        'data/chunks/chunk_{i}'
    output:
        'data/chunks/chunk_{i}.bed'
    shell:
        '''
        bedtools merge -i {input} >{output}
        '''
        

rule all_define_boostraps:
    input:
        expand('data/chunks/chunk_{i}.bed', i=DEFINED_BOOSTRAPS)
