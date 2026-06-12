#%% Extract hotspot and coldspot regions from the reference genome

import pandas as pd
import random 

from Bio import SeqIO
from Bio import SeqRecord

# Snakemake inputs
hotspots = pd.read_csv(snakemake.input.hotspots, sep='\t')
regions = pd.read_csv(snakemake.params.regions, sep='\t')
reference = SeqIO.to_dict(SeqIO.parse(snakemake.params.reference, "fasta"))
hotspot_threshold = int(snakemake.params.hotspot_threshold)
outdir = snakemake.params.outdir

# Sort hotspots according to nbirths (descending)
hotspots = hotspots.sort_values(by='nbirths', ascending=False)

# Extract sequences
hotspots_seqs = []
coldspot_seqs = []

for i, row in hotspots.iterrows():
    
    region = row['region']
    nbirths = row['nbirths']
    
    # Skip end of genome
    if region == 'mtbc0_004158;mtbc0_000001':
        continue
    
    regioninfo = regions[regions['region_name'] == region]
    regionstart = regioninfo['start'].iloc[0]
    regionend = regioninfo['end'].iloc[0]
    
    # Skip regions that are too short
    if regionend-regionstart < 10:
        continue
        
    regionseq = reference['MTBC0'].seq[regionstart:regionend]
    regionseq = SeqRecord.SeqRecord(regionseq, id=region)
    
    if nbirths >= hotspot_threshold:
        hotspots_seqs.append(regionseq)
    else:
        coldspot_seqs.append(regionseq)
        
print(f'{len(hotspots_seqs)} hotspots, {len(coldspot_seqs)} coldspots')


# Generate random genomic background sequences of same lengths as hotspots
random_seqs = []
for seq in hotspots_seqs:
    seqlen = len(seq.seq)
    seqstart = random.randint(0, len(reference['MTBC0'].seq))
    seqend = seqstart + seqlen
    random_seq = reference['MTBC0'].seq[seqstart:seqend]
    random_seq = SeqRecord.SeqRecord(random_seq, id=seq.id)
    random_seqs.append(random_seq)


# Write sequences
SeqIO.write(hotspots_seqs, f'{outdir}/hotspot_regions.fasta', "fasta")
SeqIO.write(coldspot_seqs, f'{outdir}/coldspot_regions.fasta', "fasta")
SeqIO.write(random_seqs, f'{outdir}/random_regions.fasta', "fasta")

