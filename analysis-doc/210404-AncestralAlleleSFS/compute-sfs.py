"""
Compute unfolded SFS from VCF file:
usage:
	python <VCF> <AA> <output>
Note:
	This script migth only work inside this dir
	due to relative paths used.
"""
import sys
sys.path.append("../..")
from mxbgenomes import sfs
from mxbgenomes.utils import load_populations_info

vcf_file = sys.argv[1]
aa_file = sys.argv[2]
outfile = sys.argv[3]

popinfo = load_populations_info("../../")

subpops = (
    popinfo
    .groupby('Subpopulation')
    .apply(lambda x: x.Samplename.to_list())
    .to_dict()
)

spectrum = sfs.sfs_unfolded(vcf_file, aa_file, subpops=subpops, project_haplod_size=50)
spectrum = sfs.sfs_to_frame(spectrum)
spectrum.to_csv(outfile, index=False)
