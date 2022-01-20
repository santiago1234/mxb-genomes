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


rule remove_CpGs_sites:
    input:
        vcf = 'data/jSFS/tmp/all-variants.vcf.gz',
        genome = '/data/users/smedina/data-resources/genomes/GRCh38.primary_assembly.genome-nochrprefix.fasta'
    output:
        'data/jSFS/tmp/allNoCpGs.vcf.gz',
        'data/jSFS/tmp/allNoCpGs_remove_sites.txt',
        'data/jSFS/tmp/allNoCpGs.vcf.gz.tbi'
    threads: 10
    shell:
        '''
        python ../../../scripts/removeCpGsites.py \
            {input.vcf} {input.genome} {output[0]} {output[1]}
        bcftools index --threads {threads} -t {output[0]}
        '''


rule vcf_by_boostrap_region:
    input:
        'data/jSFS/tmp/allNoCpGs.vcf.gz',
        'data/chunks/chunk_{i}.bed'
    output:
        'data/jSFS/tmp/boostraped_vcfs/all_chunk_{i}.vcf.gz'
    shell:
        '''
        bcftools view -R {input[1]} {input[0]} -Oz -o {output}
        '''


rule get_var_cat_genotypes:
    input:
        'data/jSFS/tmp/boostraped_vcfs/all_chunk_{i}.vcf.gz'
    output:
        'data/jSFS/tmp/boostraped_vcfs/all_chunk_{i}_cat_{varcat}.vcf.gz'
    params:
        v_category = vep_get_consequence_SO_term
    shell:
        '''
        python ../../../scripts/subset-vcf-by-variant-category.py \
            {input} {params.v_category} {output}
        '''


rule define_pops_to_compute_jSFS:
    output:
        'data/popsinfo.csv'
    shell:
        'python scripts/define-pops.py'


rule compute_jSFS:
    # NOTE: The jSFS is not polarized.
    input:
        vcf = 'data/jSFS/tmp/boostraped_vcfs/all_chunk_{i}_cat_{varcat}.vcf.gz',
        poplabs = 'data/popsinfo.csv'
    output:
        'data/jSFS/spectrums/spectrum_chunk_{i}_cat_{varcat}.pkl'
    shell:
        '''
        python ../../../scripts/jsfs-nonPolarized.py {input.vcf} {input.poplabs} {output}
        '''


rule compress_jSFS:
    # To save disk space.
    input:
        'data/jSFS/spectrums/spectrum_chunk_{i}_cat_{varcat}.pkl'
    output:
        'data/jSFS/spectrums/spectrum_chunk_{i}_cat_{varcat}.pkl.gz'
    shell: 'gzip {input}'


