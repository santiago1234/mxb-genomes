"""
Filter variants in VCF by VEP consequence (e.g. synonymous_variant, missense_variants, etc.)

Usage:
    python subset-vcf-by-variant-category.py <vcf_file> <category> <output_vcf>

Args:
    vcf_file: inpute VCF file
    category: SO term used by VEP to annotate variants: https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html
    out_vcf: output name for variants with the given category.

NOTES:
    input vcf should have VEP annotation.
"""

import sys
from cyvcf2 import VCF, Writer


vcf_file, category, output_vcf = sys.argv[1:]

vcf = VCF(vcf_file)
w = Writer(output_vcf, vcf)

for variant in vcf:
    vep_csq = variant.INFO.get('CSQ')
    if category in vep_csq:
        w.write_record(variant)

w.close()
vcf.close()
