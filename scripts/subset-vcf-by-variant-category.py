"""
Filter variants in VCF by VEP consequence (e.g. synonymous_variant, missense_variants, etc.)

Usage:
    python subset-vcf-by-variant-category.py <vcf_file> <category> <output_vcf>

Args:
    vcf_file: inpute VCF file
    category: SO term used by VEP to annotate variants: https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html
        multiple categries are allowed in this case use commas to separate them, example: missense_variant,synonymous_variant,intronic_variant
    out_vcf: output name for variants with the given category.

NOTES:
    input vcf should have VEP annotation.
"""

import sys
from cyvcf2 import VCF, Writer


vcf_file, category, output_vcf = sys.argv[1:]

category = category.split(',')

so_terms = ['transcript_ablation',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'start_lost',
    'transcript_amplification',
    'inframe_insertion',
    'inframe_deletion',
    'missense_variant',
    'protein_altering_variant',
    'splice_region_variant',
    'incomplete_terminal_codon_variant',
    'start_retained_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'coding_sequence_variant',
    'mature_miRNA_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'non_coding_transcript_exon_variant',
    'intron_variant',
    'NMD_transcript_variant',
    'non_coding_transcript_variant',
    'upstream_gene_variant',
    'downstream_gene_variant',
    'TFBS_ablation',
    'TFBS_amplification',
    'TF_binding_site_variant',
    'regulatory_region_ablation',
    'regulatory_region_amplification',
    'feature_elongation',
    'regulatory_region_variant',
    'feature_truncation',
    'intergenic_variant'
]


for cat in category:
    if cat not in so_terms:
        raise ValueError(f'Category {cat} not in valid so_terms')

def look_if_variant_is_cat(cat, vep_csq):
    is_cat = [x in vep_csq for x in cat]
    is_cat = sum(is_cat) > 0
    return is_cat


vcf = VCF(vcf_file)
w = Writer(output_vcf, vcf)

for variant in vcf:
    vep_csq = variant.INFO.get('CSQ')
    if look_if_variant_is_cat(category, vep_csq):
        w.write_record(variant)

w.close()
vcf.close()
