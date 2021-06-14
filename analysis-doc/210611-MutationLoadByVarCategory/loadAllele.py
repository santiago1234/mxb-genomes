import pandas as pd
import numpy as np




def l_fitness_reduction(q, h=0.5, s=0.08):
    """
    Fitnness reduction for allele
    Args:
        q: is the derived frequency of the mutation in population
        h: the dominance coefficient
        s: is the selection coefficient against the mutation
    """
    l = s * (2*h*q + (1 - 2*h)*q**2)
    return l


df_h = pd.DataFrame({'h': [0, 0.25, 0.5]})
df_s = pd.DataFrame({'s': np.arange(0, 1, 0.1)})
df_q = pd.DataFrame({'q': np.arange(0, 1, 0.01)})

d = df_h.merge(df_s, how='cross').merge(df_q, how='cross')

d['l'] = l_fitness_reduction(d.q, d.h, d.s)
d.to_csv('data/exploreload.py', index=False)
