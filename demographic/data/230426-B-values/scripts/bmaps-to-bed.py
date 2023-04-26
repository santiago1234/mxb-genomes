"""
Input data:
        
The following maps are formatted [B, length], 
after McVicker et al., 2009, where B is the 
B value (scaled from 0-1000) covering a given 
chromosome segment and length is the length 
of the segment in base pairs. Segments sum to
the length of each chromosome in the hg19 
build.

For example, the row:
 
	717 29307

assigns a B value of 0.717 (after rescaling)
to a segment of 29307 base pairs.

The output is a bed file with the following

chromosome, start, end, B value

usage:
    python bmaps-to-bed.py <infile> <outfile>

args:
    infile: bmap file
    outfile: bed file
"""

import pandas as pd
import numpy as np
import sys

infile, output = sys.argv[1:]


def bfile_to_bed(infile):
    """
    Convert a bmap file to a bed file
    """
    # ectract chromsome from filename
    chrom = infile.split('/')[-1].split('.')[0]

    bscore = pd.read_csv(infile, sep=' ', header=None)
    bscore.columns = ['B', 'length']

    # end position is the cumulative sum
    bscore['end'] = bscore['length'].cumsum()
    bscore['start'] = bscore['end'] - bscore['length']
    bscore['chrom'] = chrom

    return bscore[['chrom', 'start', 'end', 'B']]


bed = bfile_to_bed(infile)
bed.to_csv(output, sep='\t', index=False, header=False)

