#!/usr/bin/env python3
#%%
import pandas
import random
import sys

try:
    metadata = sys.argv[1]
    lineage = sys.argv[2]

except IndexError:
    sys.exit('Usage: get_lineage_strains.py <metadata> <lineage>')

md = pandas.read_csv(metadata, sep='\t')
md['LINEAGE_x'].value_counts()

# Define close outgroup for each lineage
outgroup_d = {
    'L1': 'L7',
    'L2': 'L4',
    'L3': 'L4',
    'L4': 'L2',
    'L5': 'L6',
    'L6': 'L10',
    'L7': 'L1',
    'L8': 'L5',
    'L9': 'L10',
    'L10': 'L6',
    'La1': 'La2',
    'La2': 'La1',
    'La3': 'La1',
    'C':'L8',
    'D': 'L8'
}

# Get G numbers for lineage
gnumbers = md['GNUMBER'][md['LINEAGE_x'] == lineage].to_list()

# Get random ougroup genome
outgroup = outgroup_d[lineage]
ougroup_gnumber = md['GNUMBER'][md['LINEAGE_x'] == outgroup].sample(1).to_string(index=False)

# Write to stdout
for gnumber in gnumbers:
    row = f'{gnumber}\t{lineage}\n'
    sys.stdout.write(row)
sys.stdout.write(f'{ougroup_gnumber}\t{outgroup}\n')
