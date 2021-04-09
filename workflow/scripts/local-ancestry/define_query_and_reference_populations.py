# this script is internally used in the workflow
import pandas as pd
import pickle
import numpy as np

# uncomment this line to test the script interctively
#pickle_out = open('snakemake.pickle', 'wb')
#pickle.dump(snakemake, pickle_out)
#pickle_out.close()
#snakemake = open('snakemake.pickle', 'rb')
#snakemake = pickle.load(snakemake)

# script parameters
n_pops =  snakemake.params[0]
sample_map_file =  snakemake.output[0]
query_map_file = snakemake.output[1]


#load input data

onetgp = pd.read_table("resources/1TGP-samples-meta-data/integrated_call_samples_v3.20130502.ALL.panel")
# samples from 1TGP that will be used as part of the NAT panel
onetgp_nat = np.loadtxt("resources/1TGP-samples-meta-data/native-american.txt", dtype=str)
mxb = pd.read_table("resources/genomes-metadata/50Genomes_info.txt")

# Tidy the data
mxb = mxb.loc[:, ['Sample_ID']].rename({'Sample_ID': 'Sample'}, axis=1)
mxb['Population'] = 'NAT'
mxb['Population_code'] = 'MXB'
mxb['Sample'] =  mxb.Sample.str.replace(pat='MXB', repl='MXB_')

new_names_for1tgp = {'sample': 'Sample', 'super_pop': 'Population', 'pop': 'Population_code'}
onetgp = (
        onetgp
        .rename(new_names_for1tgp, axis=1)
        .loc[:,new_names_for1tgp.values()]
    )


def native_american_panel():
    """
    Create the samples that will be the
    native american panel:
    50 MXB + onetgp_nat
    """
    natives_1tgp = onetgp[onetgp.Sample.isin(onetgp_nat)].copy()
    natives_1tgp['Population'] = 'NAT'
    return pd.concat([mxb, natives_1tgp])

def european_panel():
    """
    the european panel
    """
    sub_pops = ['IBS', 'GBR']
    return onetgp[onetgp.Population_code.isin(sub_pops)].copy()

def african_panel():
    sub_pops = ['YRI']
    return onetgp[onetgp.Population_code.isin(sub_pops)].copy()

def asian_panel():
    sub_pops = ['CHB']
    return onetgp[onetgp.Population_code.isin(sub_pops)].copy()

def query_sample_map():
    sub_pops = ['PEL', 'PUR', 'MXL', 'CLM']
    amr = onetgp[onetgp.Population_code.isin(sub_pops)].copy()
    return pd.concat([amr, mxb])


def ref_sample_map(npops):

    if str(npops) == '3':
        ref = [
                native_american_panel(),
                european_panel(),
                african_panel()
            ]
        return pd.concat(ref)

    if str(npops) == '4':
        ref = [
                native_american_panel(),
                european_panel(),
                african_panel(),
                asian_panel()
            ]
        return pd.concat(ref)

    else:
        raise ValueError('invalid number of pouplations')


ref_sample_map(npops=n_pops).to_csv(sample_map_file, index=False)
query_sample_map().to_csv(query_map_file, index=False)
