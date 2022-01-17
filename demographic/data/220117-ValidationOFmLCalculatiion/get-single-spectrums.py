import moments
import pickle
import gzip

Population = 'IBS'

sf = "../210804-Compute-jSFS/data/spectrums/5d-csq-synonymous-spectrum.pkl.gz"

with gzip.open(sf, "rb") as f:
    sf = pickle.load(f)
