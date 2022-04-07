"""
usage:
    python process-vep-output.py <inputfile> <outputfile>
VEP produces duplicated output. We do not want
to double count mutations because it will produce an
incorrect estimate.

Here, I remove duplicated elements. In case 
where a variant is assigned to multiple categories I pick
the one with the highest impact (coding only).

The idea is to have each position and variant allele represented only once.
"""

import sys
import pandas as pd
import itertools


vep_fp, vep_uniq_out = sys.argv[1:]

col_names = ['Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature_type', 'Consequence', 'Codons']
vep = pd.read_csv(vep_fp, sep='\t', comment='#', names=col_names)
print('done reading')

# Remove duplicated columns

vep = vep.drop_duplicates()


# Which variants have multiple consequences?
vep_grp = vep.groupby(['Location', 'Allele'])

vep_multiple = vep_grp.filter(lambda x: len(x) > 1)
vep_uniq = vep_grp.filter(lambda x: len(x) == 1)

vep_multiple['csq_list'] = vep_multiple.Consequence.str.split(',')
vep_uniq['csq_list'] = vep_uniq.Consequence.str.split(',')


# We want the most severe consequences (when multiple) accoring to the next rank
CSQs = [
    'stop_gained',
    'stop_lost',
    'start_lost',
    'missense_variant',
    'synonymous_variant'
]



def get_most_severq_seq(csq_list):
    """
    Args:
        csq_list: a list of predicted consequences.
    Adds the most sever from all
    the possible consequences.
    """
    #NOTE: when i run the uniq it should
    #be done based on the 3-columns
    #sometimes there are overlapping genes
    #in this case we still want only the
    #worst Consequence

    if len(csq_list) == 1:
        return csq_list[0]

    most_sever = csq_list[0]

    for x in CSQs:
        if x in csq_list:
            most_sever = x
            return most_sever
    return most_sever


vep_uniq['Consequence'] = vep_uniq.csq_list.apply(get_most_severq_seq)
vep_uniq = vep_uniq.drop(columns=['csq_list'])


vep_multiple['Consequence'] = vep_multiple.csq_list.apply(get_most_severq_seq)
vep_multiple = vep_multiple.loc[:, vep_uniq.columns]


RANK = pd.DataFrame(enumerate(CSQs), columns=['rank', 'Consequence'])

vep_m_ranked = vep_multiple.merge(RANK, how='inner', on='Consequence')

vg = vep_m_ranked.groupby(['Location', 'Allele'])


def get_best_ranked(grp):
    return grp.sort_values('rank').iloc[[0], ]

vep_m_msevere = vg.apply(get_best_ranked)
vep_m_msevere = vep_m_msevere.loc[:, vep_uniq.columns].reset_index(drop=True)


## COMBINE RESULTS
vep = pd.concat([vep_uniq, vep_m_msevere])

vep['position'] = vep.Location.str.split(':').map(lambda x: int(x[1]))

vep = vep.sort_values(by='position').drop(columns=['position'])
vep.to_csv(vep_uniq_out, sep='\t', index=False)

