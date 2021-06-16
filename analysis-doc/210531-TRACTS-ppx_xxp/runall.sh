#!bin/bash

echo "fitting tracks model"
python taino_ppx_xxp-CLM.py 1

echo "fitting tracks model"
python taino_ppx_xxp-PUR.py 1

echo "fitting tracks model"
python taino_ppx_xxp-MXL.py 1

echo "fitting tracks model"
python taino_ppx_xxp-PEL.py 1

# make dataframes for ploting
python load_results.py

