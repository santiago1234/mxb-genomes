"""
Compute the Boostraped joint-SFS.
We do this for each variant category and each chunk.
"""

rule merge_all_chroms:
    input:
        expand('../../../results/data/210713-HardyW-filters/1TGP_and_50MXB-chr{chrn}-snps-vep-mask-HW-GRCh38.vcf.gz', chrn=CHROM)
    output:
        'data/jSFS/tmp/all-variants.vcf.gz',
        'data/jSFS/tmp/all-variants.vcf.gz.tbi'
    shell:
        '''
        bcftools concat {input} -Oz -o {output}
        bcftools index -t {output}
        '''


rule vcf_by_boostrap_region:
    input:
        'data/jSFS/tmp/all-variants.vcf.gz',
        'data/chunks/chunk_{i}.bed'
    output:
        temp('data/jSFS/tmp/boostraped_vcfs/all_chunk_{i}.vcf.gz')
    shell:
        '''
        bcftools view -R {input[1]} {input[0]} -Oz -o {output}
        '''

rule remove_CpGs_sites:
    input:
        vcf = 'data/jSFS/tmp/boostraped_vcfs/all_chunk_{i}.vcf.gz',
        genome = '/data/users/smedina/data-resources/genomes/GRCh38.primary_assembly.genome-nochrprefix.fasta'
    output:
        'data/jSFS/tmp/boostraped_vcfs/allNoCpGs_chunk_{i}.vcf.gz',
        'data/jSFS/tmp/boostraped_vcfs/allNoCpGs_chunk_{i}_remove_sites.txt'
    shell:
        '''
        python ../../../scripts/removeCpGsites.py \
            {input.vcf} {input.genome} {output[0]} {output[1]}
        '''

rule get_var_cat_genotypes:
    input:
        'data/jSFS/tmp/boostraped_vcfs/allNoCpGs_chunk_{i}.vcf.gz',
    output:
        'data/jSFS/tmp/boostraped_vcfs/all_chunk_{i}_cat_{varcat}.vcf.gz'
    params:
        v_category = vep_get_consequence_SO_term
    shell:
        '''
        python ../../../scripts/subset-vcf-by-variant-category.py \
            {input} {params.v_category} {output}
        '''

