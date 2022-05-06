rule define_regions:
    output:
        'data/regions-to-sample.csv'
    shell:
        'python scripts/defined-sampled-regions.py'


rule sample_regions:
    input:
        'data/regions-to-sample.csv'
    output:
        expand('data/samples/region_region_{i}.bed', i=list(range(1, N_REGIONS + 1)))
    params:
        N = N_REGIONS
    shell:
        """
        sed '1d' {input} |\
            shuf -n {params.N} >tmp-sr.txt
        for i in {{1..{params.N}}}
        do
            head -n {params.N} tmp.txt | tail -n1 |cut -f1,2,3 >data/samples/region_region_${{i}}.bed
        done
        rm -f tmp-sr.txt
        """
