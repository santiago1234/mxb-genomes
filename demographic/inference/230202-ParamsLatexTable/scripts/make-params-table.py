"""
This script makes a table with all the parameters
in a nice format.
"""
import click

import pandas as pd


def par_cat(param):
    """
    Categorize parameter

    Args:
        param (str): parameter name
    """
    time_par = ['TA', 'TB', 'TF', 'TN']

    popsize_par = ['Ne', 'N_A', 'NB', 'NEu0',
                   'NEuF', 'NAs0', 'NAsF', 'NCmxI', 'NCmxF']

    if param in time_par:
        return 'Time event  (Thousands of years ago)'
    elif param in popsize_par:
        return 'Effective population size (# of individuals)'
    else:
        return 'Migration rate (Fraction of individuals per generation moving between populations)'


def load_infpars(ooa_pars: str, nat_pars: str, snp_cat: str):
    """Load params

    Load inferred demographic parameters

    Args:
        ooa_pars: OOA parameters, path to file
        nat_parsself: NAT parameters, path to file
        snp_cat: category from wich the give parameter where inferred
    """
    ooa_pars = pd.read_csv(ooa_pars, sep='\t')
    nat_pars = pd.read_csv(nat_pars, sep='\t')
    all_pars = pd.concat([ooa_pars, nat_pars])
    # There is a NAN value but represents the AMH
    # expansion size, here i fix this
    all_pars = all_pars.fillna('N_A').rename(columns={'#param': 'param_name'})
    all_pars['parma_type'] = all_pars.param_name.map(par_cat)
    all_pars['snp_category'] = snp_cat
    return all_pars


@click.command()
@click.option('--ooa_pars', required=True)
@click.option('--nat_pars', required=True)
@click.option('--snpcat', required=True)
@click.option('--outfile', required=True)
def main(ooa_pars, nat_pars, snpcat, outfile):
    """
    Convert vcf file to h5 format
    Args:
        ooa_pars: path to file with OOA parameters
        nat_pars: path to file with NAT parameters
        snpat: From which SNP category where
            the parameters inferred?
        outfile: File path to save resulting table
    """
    all_pars = load_infpars(ooa_pars, nat_pars, snpcat) 
    all_pars.to_csv(outfile, index=False)


if __name__ == '__main__':
    main()
