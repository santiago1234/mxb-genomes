"""
See: https://myersgroup.github.io/relate/input_data.html

Generate the recombination map needed to run relate. The input
is the usual genetic map and i compute the recombination rate
as indicated in the relate documentation.

Position (b) [integer]
Recombination rate (cM/Mb) [float]
Genetic position (cM) [float]
Denoting the ith entry of the three columns by p[i], r[i], rdist[i], the following equation holds
r[i] = (rdist[i+1] - rdist[i])/(p[i+1] - p[i]) * 1e6
pos COMBINED_rate Genetic_Map
0        2.8074 0.4103
2529     2.7778 0.4174
2601     2.9813 0.4176
This file is space delimited.

usage:
    python add_recombination_rate_to_gm.py <genetic_map> <chrn> <output>
"""
import pandas as pd
import numpy as np
import sys

gm, chrn, output = sys.argv[1:]
gm = pd.read_table(gm)


pos = gm['pos'].to_numpy()
gpos = gm['cM'].to_numpy() 
rr = (gpos[1:] - gpos[:-1]) / (pos[1:] - pos[:-1]) * 1e6

d = {
    'pos': pos[:-1],
    'COMBINED_rate': rr,
    'Genetic_Map': gpos[:-1]
}

recomb = pd.DataFrame(d)
recomb.to_csv(output, index=False, sep=' ')

