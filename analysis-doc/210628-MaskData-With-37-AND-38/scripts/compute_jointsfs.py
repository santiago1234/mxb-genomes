"""
Compute the joint SFS between populations from MXB and MXL.
usage:
    python compute_jointsfs.py <vcf_file> <aa_file> <rootToMXB> <outfile>
"""

import pandas as pd
import sys
sys.path.append("../../")
from mxbgenomes.sfs import joint_sfs
from mxbgenomes.utils import load_populations_info

# input data

vcf_file = sys.argv[1]
aa_file = sys.argv[2]
root_to_mxb = sys.argv[3]
pop_info = load_populations_info(root_to_mxb)
populations = ['MXL', 'MXB']
outfile = sys.argv[4]

spectrum, pop_index = joint_sfs(vcf_file, aa_file, populations, pop_info)


d0, d1 = spectrum.shape

cols = [pop_index[1] + '-' + str(x) for x in range(d1)]
rows = [pop_index[0] + '-' + str(x) for x in range(d0)]

spectrum = pd.DataFrame(spectrum, columns=cols, index=rows)
spectrum.reset_index().to_csv(outfile, index=False)

