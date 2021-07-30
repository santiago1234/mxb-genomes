import pickle
import numpy as np
import moments

deme_graph = "gutenkunst_2019.yaml"
options = "ooa_options.yaml"
data = "data/spectrum.pkl"
filehandler = open(data, 'rb')
data = pickle.load(filehandler)
data = data.project([40, 40, 40])

output = "./data/msl_best_fit_model.yml"
ret = moments.Demes.Inference.optimize(
    deme_graph, options, data, output=output,verbose=1, overwrite=True)



