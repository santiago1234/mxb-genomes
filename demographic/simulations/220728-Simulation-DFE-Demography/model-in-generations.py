import demes

model = '../220422-fwdpy11-initial-test/ADMIXTURE-MXL.yml'
model = demes.load(model)
model_gen = model.in_generations()

demes.dump(model_gen, 'ADMIXTURE-MXL.yml')
