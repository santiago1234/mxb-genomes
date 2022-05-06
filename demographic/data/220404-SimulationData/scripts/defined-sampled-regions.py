'''
I only sample regions where the masked sites
is at most 20%.

This script divides each chromosome in 1Mb fragments, computes
the proportion of the fragments that is masked.

I output a table that has the 1Mb regions excluding those with
more than 20% masksed sites.
'''
import pandas as pd
import pybedtools


chr_sizes = '../220113-ConstructBoostrapedDatasets/data/human-autosomes.genome'

chr_sizes = pd.read_csv(chr_sizes, sep='\t')

def make_1mb_windows_in_chr(chro):
    chromosome = chr_sizes[chr_sizes.chrom == chro]
    chrsize = chromosome.iloc[0, 1]
    region_size = int(1e6)
    current_pos = 0
    regions = []

    chromname = 'chr' + str(chro)
    while True:
        if (current_pos + region_size) < chrsize:
            regions.append((chromname, current_pos, current_pos + region_size))
            current_pos += region_size
        else:
            break
    
    return pd.DataFrame(regions, columns=['chrn', 'start', 'end'])


mask = pybedtools.BedTool('../../../results/data/210305-merged-with-1TGP/strict-mask/20160622.allChr.mask.bed')

def region_percent_mask(chrn, start, end):
    print(f'Getting {chrn} {start}, {end}')
    br = pybedtools.BedTool(f'{chrn} {start} {end}', from_string=True)
    br_mask = br.intersect(mask)

    lens = [x.end - x.start for x in br_mask]
    return 1 - sum(lens) / (end - start)
        

## TODO: Modigfy the number below to run all genome
regions = pd.concat([make_1mb_windows_in_chr(i) for i in range(21, 23)])

## add a column with the percentage masked
regions['masked_percent'] = regions.apply(lambda x: 
        region_percent_mask(*(x.chrn, x.start, x.end)), axis=1)

## remove regions with more than 80% mask
regions = regions[regions.masked_percent < 0.20]

## save regions
regions.to_csv('data/regions-to-sample.csv', sep='\t', index=False)
