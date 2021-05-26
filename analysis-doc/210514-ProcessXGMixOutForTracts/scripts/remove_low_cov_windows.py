"""
Drop windows with low SNP density
usage:
    python remove_low_cov_windows.py <input> <N> <output>
Args:
    input: The input msp file containing the ancestry assigments.
    N: The number of windows to drop, the top <N> windows with
        lowest SNP density will be dropped.
    output: file name for output file
"""
import pandas as pd
import sys


# input parameters
msp_file = sys.argv[1]
drop_topN = int(sys.argv[2])
out_file = sys.argv[3]

msp = pd.read_csv(msp_file)

# compute SNP density

window_length = msp.epos - msp.spos
msp['snp_density'] = msp['n snps'] / window_length

drop_low_density = (
    msp
    .sort_values('snp_density')  # sort by snp density
    .iloc[drop_topN:]
    .sort_values(['#chm', 'spos'])  # re-sort by position
    .drop('snp_density', axis=1)
)

drop_low_density.to_csv(out_file, index=False)
