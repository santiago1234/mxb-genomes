# Merge 50 genomes with 1TGP high coverage data


rule subset_pop_list:
    # this rule generates a list with the population to
    # subset from 1TGP data
    input:
        pop_meta = "resources/1TGP-samples-meta-data/igsr-1000genomes.tsv"
    output:
        "results/data/210305-merged-with-1TGP/pops-to-subset.txt"
    params:
        pops = '(CHB|YRI|IBS|BGR|MXL|PEL|CLM|PUR)' # use a regular expression
    shell:
        """
        cut -f1,4 {input} |\
            grep -E '{params.pops}' |\
            cut -f1 >{output}
        """

    
