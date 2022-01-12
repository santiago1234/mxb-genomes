import pandas as pd

# I will use chr22 

variants = '../data/tmp/variants/chr22.bed.gz'
variants = pd.read_csv(variants, sep='\t', names=['chr', 'spos', 'epos', 'var', 'strand', 'id'])

counts = variants.spos.value_counts().value_counts() 
print(counts)
