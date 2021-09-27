# generate the vcf with the synonymous variants
snakemake -j1 data/vcf-synonymous.vcf.gz
snakemake -j1 data/all-ancestral-allel.csv
