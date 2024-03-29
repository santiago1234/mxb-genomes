"""
Take the genome divide in non-overlapping windows
and split in files, such that each file spans a max_length
of 7.5 Mb.

The expected length of each bootstrap chunk will be ~6Mb.
This is considering that some regions of the genome are maked.
"""
import pandas as pd

d = pd.read_csv('data/chunks/genome-in-1MB-non-overlaping-windows.bed', sep='\t', names=['chr', 'spos', 'epos'])
outpath = 'data/chunks/'
d['width'] = d.epos - d.spos

max_length = 10000000
length_chunk = 0
chunks = list()
chunks.append([])

i = 0
for index, row in d.iterrows():
    length_chunk += row.width
    
    if length_chunk < max_length:
        chunks[i].append(index)
    
    else:
        length_chunk = 0
        chunks.append([])
        i += 1


# Now save each boostrap chunk

N_chunks = len(chunks)

for i in range(N_chunks):
    print(f'Saving chunk {i+1} ...')
    outfile = f'data/chunks/chunk_{i+1}'
    d_chunk = d.iloc[chunks[i], [0, 1, 2]]
    d_chunk.to_csv(outfile, sep='\t', index=False, header=False)

