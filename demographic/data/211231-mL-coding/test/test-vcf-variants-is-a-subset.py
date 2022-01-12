import pandas as pd

universe_missense = '../data/tmp/varcat/chr22-missense.txt'
universe = pd.read_csv(universe_missense, sep='\t', names=['main', 'pos', 'alt', 'gene', 'codon'])

genotype_vars = 'tmp/genotypes-variants.txt'
genotype_vars = pd.read_csv(genotype_vars, sep='\t', names=['chrn', 'pos', 'ref', 'alt'])

# make the names to match
universe['ref'] = universe.main.str[-3]
universe['pos'] = universe.pos.str.split(':').apply(lambda x: int(x[-1]))


# The universe positons
pos_universe = set(universe.pos.values)

positions_not_found = []
for index, row in genotype_vars.iterrows():
    if not (row.pos in pos_universe):
        print(f'Test failed at position: {row.pos}')
        positions_not_found.append(row.pos)

f = open('regions-not-found.bed', 'w')

for p in positions_not_found:
    f.write(f'22\t{p - 1}\t{p + 1}\n')
f.close()

