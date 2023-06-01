"""
Combine stats per region into a single file.
"""
import sys
import pickle
from glob import glob

POPS = ['YRI', 'MXL', 'MXB', 'IBS', 'CHB']


def load_ld_stat_pop(population):
    # termination is txt but it's actually a pickle
    files = glob(f'results/ld_stats/{population}*.txt')

    # get chunk id from file name: ld_stats/MXB-region69-ld_stats.txt
    def get_chunk_id(file_name):
        reg = file_name.split('/')[-1].split('-')[1]
        return int(reg.replace('region', ''))


    data = {}

    for file in files:
        chunk_id = get_chunk_id(file)
        with open(file, 'rb') as f:
            data[chunk_id] = pickle.load(f)
    
    return data

data = {}
for pop in POPS:
    data[pop] = load_ld_stat_pop(pop)

with open('results/all-pop-ld-stats.pkl', 'wb') as f:
    pickle.dump(data, f)
