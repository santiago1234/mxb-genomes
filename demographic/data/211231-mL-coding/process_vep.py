import pandas as pd

csq = 'missense_variant'
outname = csq + '.txt'
cols = "Uploaded_variation,Location,Allele,Gene,Feature_type,Consequence,Codons".split(',')
vep = pd.read_csv('data/variant_effect_output.txt', sep='\t', comment='#', names=cols)

vep = vep[vep.Consequence.str.contains(csq)]
cols.remove('Consequence')

vep = vep.drop_duplicates(subset=cols)

vep['context'] = vep.Uploaded_variation.str.extract(r'_([ACTG]*)_')
vep['context'] = vep.context.str[1:-1]

vep.to_csv(outname, sep='\t', index=False)
