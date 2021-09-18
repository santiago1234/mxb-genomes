import pandas as pd
import glob

results = [pd.read_csv(x) for x in glob.glob('data/counts/derived-counts-individual*')]
results = pd.concat(results, ignore_index=True)

gp_dev_counts = results.groupby(['Samplename', 'variant', 'Region', 'VarFreq'])['derived_count']
counts = gp_dev_counts.sum().reset_index()
counts.to_csv('results/derived_counts.csv', index=False)
