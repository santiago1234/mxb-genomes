import pybedtools
import pandas as pd

chr_sizes = '../220113-ConstructBoostrapedDatasets/data/human-autosomes.genome'

chr_sizes = pd.read_csv(chr_sizes, sep='\t')

def make_1mb_windows_in_chr(chro):
    chromosome = chr_sizes[chr_sizes.chrom == chro]
    chrsize = chromosome.iloc[0, 1]
    region_size = int(1e6)
    current_pos = 0
    regions = []

    chromname = str(chro)
    while True:
        if (current_pos + region_size) < chrsize:
            regions.append((chromname, current_pos, current_pos + region_size))
            current_pos += region_size
        else:
            break
    
    return pd.DataFrame(regions, columns=['chrn', 'start', 'end'])


exons = pybedtools.BedTool('data/exons-sorted-merged.bed')

def region_exon_len(chrn, start, end):
    print(f'Getting {chrn} {start}, {end}')
    br = pybedtools.BedTool(f'{chrn} {start} {end}', from_string=True)
    br_exons = br.intersect(exons)

    lens = [x.end - x.start for x in br_exons]
    return sum(lens)

regions = pd.concat([make_1mb_windows_in_chr(i) for i in range(1, 23)])

## add a column with the exon length
regions['exon_len'] = regions.apply(lambda x: 
        region_exon_len(*(x.chrn, x.start, x.end)), axis=1)

regions = regions.sort_values(by='exon_len', ascending=False)

# Keep the top 3rd gene rich regions 
# The genome is 3,000 Mb so 1000 regions is 1/3
regions.head(1000).to_csv('data/regions-to-sample.csv', sep='\t', index=False)
