# Define start parameters for each model and population

POPS = ['MXL', 'CLM', 'PUR', 'PEL']

# -------------------------------------------------------------------------
# This start parameters work for all populations
startparams_ppx_xxp = [
    0.0667324402332306,  # t1
    0.055764797329902666  # t2
]


STP_ppx_xxp = {x: startparams_ppx_xxp for x in POPS}

# -------------------------------------------------------------------------
startparams_ppx_xxp_pxx = [
    0.16347885,
    0.08063866,
    0.0799999,
    0.06329345
]


STP_ppx_xxp_ppx = {x: startparams_ppx_xxp_pxx for x in POPS}

# -------------------------------------------------------------------------

startparams_ccx_xxp = [
    0.0533,
    0.0533,
    0.1536,
    0.0913,
    0.11673
]


startparams_ccx_xxp_PUR = [
    0.007,
    0.0353,
    0.154,
    0.127,
    0.099
]


STP_ccx_xxp = {x: startparams_ccx_xxp for x in POPS}
STP_ccx_xxp['PUR'] = startparams_ccx_xxp_PUR

# -------------------------------------------------------------------------

# This works for MXL and PUR
startparams_ppx_ccx_xxp = [
    0.2140,  # initEu_frac
    0.0533,  # frac1
    0.0533,  # frac2
    0.1536,  # t1
    0.0913,  # frac3
    0.11673  # t2
]


startparams_ppx_ccx_xxp_PEL = [
    0.263,
    0.119,
    0.0155,
    0.144,
    0.0482,
    0.069
]


startparams_ppx_ccx_xxp_CLM = [
    0.563,
    0.011,
    0.0622,
    0.13,
    0.118,
    0.09
]


startparams_ppx_ccx_xxp_PUR = [
    0.186,
    0.005,
    0.0346,
    0.148,
    0.185,
    0.120
]


STP_ppx_ccx_xxp = {
    'MXL': startparams_ppx_ccx_xxp,
    'PUR': startparams_ppx_ccx_xxp_PUR,
    'CLM': startparams_ppx_ccx_xxp_CLM,
    'PEL': startparams_ppx_ccx_xxp_PEL
}

START_PARAMS = {
    'ppx_xxp': STP_ppx_xxp,
    'ppx_xxp_pxx': STP_ppx_xxp_ppx,
    'ccx_xxp': STP_ccx_xxp,
    'ppx_ccx_xxp': STP_ppx_ccx_xxp
    }
