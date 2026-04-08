#%% Extract hotspot and coldspot regions from the reference genome

import pandas as pd
from Bio import SeqIO


#%%
hotspots = pd.read_csv('../results/hotspot_summary.tsv', sep='\t')
regions = pd.read_csv('../../data/MTBC0/MTBC.genomic_regions.extended.tsv', sep='\t')
reference = SeqIO.to_dict(SeqIO.parse('../../data/MTBC0/MTBC0_v1.1.fasta', "fasta"))

hotspot_threshold = 10

#%%
# Sort hotspots according to nbirths (descending)
hotspots = hotspots.sort_values(by='nbirths', ascending=False)

hotspots_seqs = []
coldspot_seqs = []

for i, region in hotspots.iterrows():
    
    #  Get matching region in regions
    region = regions[regions['region'] == region['region']]
    region = region.iloc[0]
    
    
    
    if region['nbirths'] > hotspot_threshold:
            


