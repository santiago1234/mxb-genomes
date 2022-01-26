'''
Get the end time for the OOA population,
also called CHB start time.
usage:
    python get-CHB-start-time.py ooa.yml
'''
import yaml
import sys

demes_graph = sys.argv[1]
demes = yaml.safe_load(open(demes_graph, 'r'))

#Get the deme for the OOA population
ooa_deme = [d for d in demes['demes'] if d['name'] == 'OOA']
ooa_end_time = ooa_deme[0]['epochs'][0]['end_time']

print(f'  upper_bound: {ooa_end_time}')
