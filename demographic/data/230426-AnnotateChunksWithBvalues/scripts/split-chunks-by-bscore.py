"""
Split chunks by bscore, I do four categories each category
is a quantile of bscore.
"""
import pandas as pd

bscores = pd.read_csv('results/bscores.txt', sep=' ', header=None)
bscores.columns = ['chunk', 'bscore']

quartiles = bscores['bscore'].quantile([0.25, 0.5, 0.75])

def assign_quartile(b_score):
    if b_score <= quartiles[0.25]:
        return 1
    elif b_score <= quartiles[0.5]:
        return 2
    elif b_score <= quartiles[0.75]:
        return 3
    else:
        return 4

bscores['quartile'] = bscores['bscore'].apply(assign_quartile)

# Now save a file for each quartile

for i in range(1, 5):
    qb = bscores[bscores['quartile'] == i]
    qb[['chunk']].to_csv(f'data/chunks-by-bscores-q{i}.txt',
                                             sep=' ', header=None, index=None)
