"""
Compute confidence intervals using the Godambe method.

usage:
    python Estimate-CI.py <SPECTRUM_GENOME> <ML_GENOME> <BOOTSTRAPS> <PROJECT_SIZE> <MARGINALIZE_NAT> <BEST_GUEST_model> <PARAMETERS> <OUTPUT_TABLE>


Args:
    SPECTRUM_GENOME: File to the same data SFS used in the origional optimization.
    ML_GENOME: File to the sequence length-scaled mutation rate.
    BOOTSTRAPS: Bootstraped replicates.
    PROJECT_SIZE: N, interger the size to project. I recommend to 30.
    MARGINALIZE_NAT: yes or no. If yes we marginalize the MXB/NAT population.
    BEST_GUEST_model: The file path to the optimized demes graph.
    PARAMETERS: The same options file used in the original optimization.
    OUTPUT_TABLE: Path to table with uncerts.
"""
import moments
import pickle
import gzip
import sys

SPECTRUM_GENOME, ML_GENOME, BOOTSTRAPS, PROJECT_SIZE, MARGINALIZE_NAT, BEST_GUEST_model, PARAMETERS, OUTPUT_TABLE, = sys.argv[1:]


if MARGINALIZE_NAT.lower() == 'yes':
    MARGINALIZE_NAT = True
else:
    MARGINALIZE_NAT = False

PROJECT_SIZE = int(PROJECT_SIZE)

def read_spectrum_file(spec_file):
    '''Read moments.Spectrumb'''
    with gzip.open(spec_file, "rb") as f:
        sf = pickle.load(f)

    return sf


def read_mL_file(mL_file):
    '''Read mL from file'''
    f = open(mL_file, 'r')
    mL = f.readlines()[0]
    mL = mL.replace('mL:', '')
    mL = mL.strip()
    mL = float(mL)
    f.close()
    return mL


def process_spectrum(sf, projection_size, marginalize_NAT):
	"""
    Args:
        sf: moments.Spectrum
        projection_size: Size to project
	project, marginalizes, and folds spectrum
    """
	if marginalize_NAT:
		project = [projection_size] * 3
		sf= sf.marginalize([sf.pop_ids.index('MXB')])
		sf = sf.project(project)
	else:
		project = [projection_size] * 4
		sf = sf.project(project)

	return sf.fold()



data_genome = read_spectrum_file(SPECTRUM_GENOME)
mL_genome = read_mL_file(ML_GENOME)

print('loading data ...')
with gzip.open(BOOTSTRAPS, 'rb') as f:
        BOOTSTRAPS = pickle.load(f)


bootstraps = [x[1] for x in BOOTSTRAPS]
bootstraps_mL  = [x[2] for x in BOOTSTRAPS]
## apply to data
data_genome = process_spectrum(data_genome, PROJECT_SIZE, MARGINALIZE_NAT)
bootstraps = [process_spectrum(x, PROJECT_SIZE, MARGINALIZE_NAT) for x in bootstraps]

std_err = moments.Demes.Inference.uncerts(
    BEST_GUEST_model,
    PARAMETERS,
    data_genome,
    bootstraps=bootstraps,
    uL=mL_genome,
    bootstraps_uL=bootstraps_mL,
    method="GIM",
	verbose=50,
	output='results.tab'	
)

