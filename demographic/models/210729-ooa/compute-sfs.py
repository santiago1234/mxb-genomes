import numpy as np
import sys
import pickle
import moments
sys.path.append('../../../')
from mxbgenomes.sfs import joint_sfs
from mxbgenomes.utils import load_populations_info

vcf_file = "../../../results/data/210713-HardyW-filters/1TGP_and_50MXB-chr22-snps-vep-mask-HW-GRCh38.vcf.gz"
aa_file = "../../../analysis-doc/210506-AncestralAlleleData/data/aa-chr22.csv"
root_to_mxb = '../../../'
pop_info = load_populations_info(root_to_mxb)

populations = ['YRI', 'IBS', 'CHB']

spectrum, pop_index = joint_sfs(vcf_file, aa_file, populations, pop_info)

spectrum = moments.Spectrum(spectrum, pop_ids=[pop_index[i] for i in range(3)], data_folded=False)
spec_file = open('data/spectrum.pkl', 'wb')
pickle.dump(spectrum, spec_file)
spec_file.close()
