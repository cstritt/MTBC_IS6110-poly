#%%
from Bio import SeqIO
from BCBio import GFF
import os
import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Extract CDS and intergenic regions from reference genome')
    parser.add_argument('--input_fasta', required=True, help='Input reference genome FASTA')
    parser.add_argument('--input_gff', required=True, help='Input GFF annotation file')
    parser.add_argument('--output_fasta', required=True, help='Output CDS+intergenic FASTA')
    parser.add_argument('--output_burden', required=True, help='Output burden TSV file')
    return parser.parse_args()

args = parse_args()

input_fasta = args.input_fasta
input_gff = args.input_gff
output_file = args.output_fasta
output_burden_file = args.output_burden

"""
Extract CDS regions and intergenic regions from a reference genome based on GFF annotations.

Args:
    input_fasta (str): Path to the input reference genome in FASTA format
    input_gff (str): Path to the input GFF annotation file
    output_file (str): Path to the output FASTA file containing CDS and intergenic regions
"""
#%%
# First, parse the GFF file to extract CDS features and organize them by scaffold
cds_features = {}

feature_types = set()

with open(input_gff) as in_handle:
    for rec in GFF.parse(in_handle):
        scaffold_id = rec.id
        if scaffold_id not in cds_features:
            cds_features[scaffold_id] = []
            
        for feature in rec.features:
            
            feature_types.add(feature.type)
            feature_id = feature.qualifiers['ID'][0]
            
            descr = feature.qualifiers['gene'][0] if 'gene' in feature.qualifiers else ''
            
            # Store CDS features with their coordinates
            cds_features[scaffold_id].append({
                'start': min(feature.location.start, feature.location.end),
                'end': max(feature.location.start, feature.location.end),
                'strand': feature.location.strand,
                'id': feature.id if feature.id else f"{scaffold_id}_CDS_{len(cds_features[scaffold_id])+1}",
                'type': feature.type,
                'descr' : descr
            })

# Sort CDS features by their start position
for scaffold in cds_features:
    cds_features[scaffold].sort(key=lambda x: x['start'])

#%%
# Load the reference genome
genome = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))

#%%
# Prepare to write output sequences
output_records = []
burden_regions = []  # feature id\tchrom:start-end

# Process each scaffold in the genome
for scaffold_id, scaffold_seq in genome.items():
    if scaffold_id not in cds_features:
        # If no CDS annotations for this scaffold, add the entire sequence as intergenic
        output_records.append(scaffold_seq)
        continue
        
    cds_list = cds_features[scaffold_id]
    seq_length = len(scaffold_seq)
    
    # Add intergenic region before the first CDS
    first_cds = cds_list[0]
    if first_cds['start'] > 0:
        inter_seq = scaffold_seq[0:first_cds['start']]
        inter_seq.id = f";{first_cds['id'][5:]}"
        inter_seq.description = ''
        output_records.append(inter_seq)
        burden_regions.append([inter_seq.id, f"{scaffold_id}:1-{first_cds['start']}"])
    
    # Add CDS regions
    for cds in cds_list:
        cds_id = cds['id'].replace('gene-', '')
        cds_seq = scaffold_seq[cds['start']:cds['end']]
        if cds['strand'] == -1:
            cds_seq = cds_seq.reverse_complement()
        
        if cds['type'] == 'gene':
            cds_seq.id = cds_id
            cds_seq.description = cds['descr']
        else: 
            cds_seq.id = cds_id
            cds_seq.description = cds['type']
        
        output_records.append(cds_seq)
        burden_regions.append([cds_seq.id, f"{scaffold_id}:{cds['start']+1}-{cds['end']}"])
        
        # Add intergenic region between this CDS and the next
        if cds != cds_list[-1]:  # If not the last CDS
            next_cds = cds_list[cds_list.index(cds) + 1]
            cds_id = cds['id'].replace('gene-', '')
            next_cds_id = next_cds['id'].replace('gene-', '')
            
            if cds['end'] < next_cds['start']:
                inter_seq = scaffold_seq[cds['end']:next_cds['start']]
                inter_seq.id = f"{cds_id};{next_cds_id}"
                inter_seq.description = ''
                output_records.append(inter_seq)
                burden_regions.append([inter_seq.id, f"{scaffold_id}:{cds['end']+1}-{next_cds['start']}"])
    
    # Add intergenic region after the last CDS
    last_cds = cds_list[-1]
    if last_cds['end'] < seq_length:
        inter_seq = scaffold_seq[last_cds['end']:seq_length]
        inter_seq.id = f"{last_cds['id'][5:]};"
        inter_seq.description = 'intergenic'
        output_records.append(inter_seq)
        burden_regions.append([inter_seq.id, f"{scaffold_id}:{last_cds['end']+1}-{seq_length}"])

#%%
# Write the output file
with open(output_file, "w") as out_handle:
    SeqIO.write(output_records, out_handle, "fasta")

# %% Writen region file for gene burden analyis (name\tchrom:start-end)
with open(output_burden_file, "w") as out_handle:
    for region in burden_regions:
        out_handle.write(f"{region[0]}\t{region[1]}\n")
    

# %%
