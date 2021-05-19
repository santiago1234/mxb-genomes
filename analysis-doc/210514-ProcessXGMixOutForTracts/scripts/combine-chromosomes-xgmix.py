"""
Merge the XGMix output files
per chromosome into a single one
"""

import sys
sys.path.append("../../")
from mxbgenomes.localancestry import processxgmixout

path_to_xgmix_out = sys.argv[1]
outfile = sys.argv[2]

predictions = processxgmixout.load_xgmix_output(path_to_xgmix_out)
predictions.to_csv(outfile, index=False)
