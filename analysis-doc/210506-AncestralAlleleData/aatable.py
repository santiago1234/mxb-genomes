"""
Script to make a table containing the variand id
and the ancestral allele
usage:
	python <vcf_file> <outfile.csv>
"""
import pandas as pd
import allel
import sys

vcf_file = sys.argv[1]
out_file = sys.argv[2]

callset = allel.read_vcf(vcf_file)


# Alternative:

# Here is how I can extract the ancestral alle
csq = allel.read_vcf(vcf_file, fields='CSQ')


# which variants do the reference alle matches the ancestral  allele

# Make table with the following information:
# Variant_id, Ancestral allel

ancestral_allel = csq['variants/CSQ']
variant_id = callset['variants/ID']

df_ancestral = pd.DataFrame({'ID': variant_id, 'AA': ancestral_allel})

df_ancestral['AA'] = df_ancestral.AA.str.upper()
df_ancestral.to_csv(out_file, index=False)

