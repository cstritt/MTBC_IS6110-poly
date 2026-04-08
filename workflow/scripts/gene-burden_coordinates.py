#!/usr/bin/env python


#%% 

import re

go_terms = ['GO:0006310', 'GO:0006281', 'GO:0006260']
mtbc0_annot = '/scicore/home/gagneux/stritt0001/TB/projects/MTBC_IS6110-poly/data/MTBC0/MTBC0v1.1_PGAP_annot.gff'
h37rv_annot = '/scicore/home/gagneux/stritt0001/TB/projects/MTBC_IS6110-poly/data/H37Rv/GCF_000195955.2_ASM19595v2_genomic.gff'
annot_map = '/scicore/home/gagneux/stritt0001/TB/projects/MTBC_IS6110-poly/data/MTBC0/mtbc0_to_h37rv.genemap.tsv'


#%% Grep GO terms in MTBC0 annotation
r3_genes = []
with open(mtbc0_annot) as f:
    for line in f:
        for term in go_terms:
            if term in line and line not in r3_genes:
                r3_genes.append(line)
    
print(f'{len(r3_genes)} genes matching GO terms')


#%% Get corresponding orthologs in H37Rv
genemap = {}

with open(annot_map) as f:
    next(f)
    for line in f:
        fields = line.strip().split()
        genemap[fields[0]] = fields[1]
 
genemap

#%%
pattern = r'ID=([^;]+)'

rv_orthologs = []

for line in r3_genes:
    match = re.search(pattern, line)
    gene_id = match.group(1)[4:]
    
    if gene_id in genemap:
        rv = genemap[gene_id]
        rv_orthologs.append(rv)
    else:
        print(gene_id)
        
rv_orthologs
    
#%% Go through H37Rv annotation and write to files: a burden input file with R3 genes only, one with all genes
"""
gene1    chrom:start-end

"""

pattern = r'locus_tag=([^;]+)'

chromosome = 'MTB_anc'

with open(h37rv_annot) as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if fields[2] == 'gene':
            match = re.search(pattern, fields[-1])
            locus_tag = match.group(1)
            
            start = fields[3]
            end = fields[4]
            
            if locus_tag in rv_orthologs:
                print(fields)
            
        


# %%
