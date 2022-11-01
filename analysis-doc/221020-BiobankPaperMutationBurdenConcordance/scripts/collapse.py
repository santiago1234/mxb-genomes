import sys
import re
import glob

import pandas as pd

sys.path.append('../../')

from mxbgenomes.utils import load_populations_info


dfiles = glob.glob('data/counts/*csv')
output = 'results/alt-counts.csv'

def decode_meta_from_filename(afile):
    '''
    Extract the metada from the file name
    '''
    chrom = re.findall('chrn\d\d?', afile)[0]
    vcat = re.findall('cat([A-Z]*)', afile)[0]
    region = re.findall('reg([A-Z]*)', afile)[0]
    return chrom, vcat, region


def load_data(afile):
    '''
    Load counts as dataframe 
    and adds the metadata
    '''
    data = pd.read_csv(afile)
    chrom, vcat, region = decode_meta_from_filename(afile)
    data['chrom'] = chrom
    data['vcat'] = vcat
    data['region'] = region

    return data


all_data = pd.concat([load_data(x) for x in dfiles])

all_data = (
        all_data
        .groupby(['sample', 'vcat', 'region'])
        .alt_counts
        .sum()
        .reset_index()
        .rename(columns={'sample': 'Samplename'})
    )

popinfo = load_populations_info('../../')
popinfo = popinfo.loc[:, ['Samplename', 'Subpopulation']]
all_data = pd.merge(all_data, popinfo, on=['Samplename'])

all_data.to_csv(output, index=False)
