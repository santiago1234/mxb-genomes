"""
Generate the counts table. Count the number of mutations (Contex)
Output file is meant to use for computing mL.

usage: python count-context-seqs.py <input-variants> <output-file>
Args:
    <input-variants>: input file with variants.
    <output-file>: file-path
"""
import sys
import pandas as pd

infile, outputfile = sys.argv[1:]

cols = "Uploaded_variation,Location,Allele,Gene,Codons".split(',')
d = pd.read_csv(infile, sep='\t', names=cols)


# Get the context sequence
d['context'] = d.Uploaded_variation.str.extract(r'_([ACGT]{5})_')
# we want the triplet in the center
d['context'] = d.context.str[1:-1]

counts = d.groupby(['context', 'Allele']).size().reset_index()
counts['ref'] = counts.context.str[1:-1]

counts = counts.rename(columns={0: 'N', 'Allele': 'alt'}).loc[:,['context', 'ref', 'alt', 'N']]
counts.to_csv(outputfile, index=False, sep='\t')
