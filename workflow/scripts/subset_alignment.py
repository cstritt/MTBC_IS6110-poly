#!/usr/bin/env python3

#%%
import sys

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

def subset_alignment(alignment_file, sequence_names, output_file):
    # Read the alignment file
    alignment = AlignIO.read(alignment_file, 'fasta')
    
    # Read the list of sequence names
    samples = []
    with open(sequence_names, 'r') as f:
        for line in f:
            fields = line.split('\t')
            samples.append(fields[0])

    # Subset the alignment using the list of sequence names
    subset_alignment = MultipleSeqAlignment([seq for seq in alignment if seq.id in samples])

    # Remove invariable sites
    variable_sites = [i for i in range(subset_alignment.get_alignment_length()) if len(set([seq[i] for seq in subset_alignment])) > 1]
    subset_alignment = subset_alignment[:, variable_sites]

    # Write the subset alignment to a new file
    AlignIO.write(subset_alignment, output_file, 'fasta')

#%%

def main():
    
    try:
        alignment_file = sys.argv[1]
        sequence_names = sys.argv[2]
        output_file = sys.argv[3]
    except IndexError:
        sys.exit('Usage: subset_alignment.py <alignment_file> <sequence_names.txt> <output_file>')
    
    subset_alignment(alignment_file, sequence_names, output_file)
    
if __name__ == '__main__':
    main()
# %%
