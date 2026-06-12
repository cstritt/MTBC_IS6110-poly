"""

"""

import random
import sys
from pathlib import Path
import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

# ── snakemake objects ────────────────────────────────────────────────────────

min_length = int(snakemake.params.min_length)
max_length = int(snakemake.params.max_length)
max_weights = int(snakemake.params.max_weights)
max_total_seqs = int(snakemake.params.max_total_seqs)

log_path     = str(snakemake.log[0])
seed         = int(snakemake.params.seed)

random.seed(seed)

# Redirect stderr to log
import io
log_fh = open(log_path, "w")

def log(msg):
    log_fh.write(msg + "\n")
    log_fh.flush()


# ── helpers ──────────────────────────────────────────────────────────────────

def load_reference(fasta_path):
    rec = next(SeqIO.parse(fasta_path, "fasta"))
    seq = str(rec.seq).upper()
    log(f"[ref] {rec.id}  length={len(seq):,} bp")
    return seq


def load_births(tsv_path):
    """Return dict {site_id: float(births)}."""
    df = pd.read_csv(tsv_path, sep="\t")
    # tolerate either 'site_id' or first column as ID
    id_col = "site_id" if "site_id" in df.columns else df.columns[0]
    births_col = "births" if "births" in df.columns else df.columns[1]
    df[id_col] = [x[1:] for x in df[id_col].astype(str)]
    result = dict(zip(df[id_col].astype(str), df[births_col].astype(float)))
    log(f"[births] {tsv_path}  →  {len(result):,} entries")
    return result


def load_flanking_sequences(anchors_tsv):
    """
    Load ALL_anchors.tsv and store cFS.   
    Returns a dictionary[(strain, anchor_id)] = SeqRecord().
    """
    records = {}
    
    anchor_df = pd.read_csv(anchors_tsv, sep="\t", comment="#")
    
    for i, row in anchor_df.iterrows():
        anchor_id = row['anchor_id']
        strain = row['strain']
        sequence = row['consensus']
        records[(strain, anchor_id)] = SeqRecord(Seq(sequence), id=f'{strain}.{anchor_id}', name="", description="")
    
    log(f"[anchors] loaded {len(records):,} anchors from FASTA")
    return records


def join_flanking_sequences(refernce_insertions, flanking_sequences, anchor_map, 
                            min_length = 50, max_length = 150):
    """
    For each insertion with identifiable reference position and a 5'-'3 cFS pair, 
    join the 5' and 3' flanking sequences and remove the target site duplication. 
    
    Cut the flanking sequences to equal sizes, that is, to the size of the smaller. 
    If the resulting joined sequence is longer than max_length, cut it to max_length.
    
    Randomly subsample the joint_seqs to max_total_seqs, to avoid issues with STREME.
    
    
    Returns a dictionary[(strain, anchor_id)] = SeqRecord().
    In the description, add the ID of the unique flanking sequence. 
    
    Args:
        refernce_insertions (_type_): _description_
        flanking_sequences (_type_): _description_
    """
    refins = pd.read_csv(refernce_insertions, sep="\t")
    anchormap = pd.read_csv(anchor_map, sep="\t", names=['strain', 'cFS_id', 'uFS_id', 'side'])
    
    joint_seqs = defaultdict(list)  # Avoid checking if key exists
    
    anchor_lookup = {}
    for _, row in anchormap.iterrows():
        key = (row['strain'], row['cFS_id'])
        anchor_lookup[key] = row['uFS_id']
    
    for i, row in refins.iterrows():
        strain = row['strain']
        anchor_5 = row['anchor_5']
        anchor_3 = row['anchor_3']
        tsd = row['TSD']
        
        if anchor_5.startswith('3'):
            flanking_5 = anchor_3
            flanking_3 = anchor_5
        else:
            flanking_5 = anchor_5
            flanking_3 = anchor_3
            
        seq5 = flanking_sequences[(strain, flanking_5)]
        seq3 = flanking_sequences[(strain, flanking_3)]
        
        if seq5 is None or seq3 is None:
            continue
        
        # Get the name of the corresponding unique flanking sequence
        # Structure of anchor  map: strain anchor_id anchor_id_unique
        uFS_id = anchor_lookup.get((strain, flanking_5)) or anchor_lookup.get((strain, flanking_3))
        
        if uFS_id is None:
            log(f"[anchors] no unique flanking sequence found for {strain} {flanking_5} {flanking_3}")
            continue

        # Remove TSD from 5' flanking sequence
        seq5.seq = seq5.seq[:len(seq5)-len(tsd)]
        
        # Join 5' and 3' flanking sequences. Cut from beginning of 5' and end of 3'
        cut_to = min(len(seq5.seq), len(seq3.seq))
        
        # Skip small sequences
        if 2*cut_to < min_length:
            continue
        
        # Cut to max_length
        if (2*cut_to) > max_length:
            cut_to = max_length // 2

        jointseq = SeqRecord(
            seq5.seq[-cut_to:] + seq3.seq[:cut_to], 
            id=f'{uFS_id}_{strain}', 
            name="", 
            description=f'{flanking_5}_{flanking_3}'
            )
        
        joint_seqs[uFS_id].append(jointseq)
    
    log(f"[anchors] joined {len(joint_seqs):,} flanking sequences")
    return joint_seqs
    

def write_motif_fasta(joint_seqs, weights, 
                      max_weight = 50, max_total_seqs = 5000):
    """ 
    Write FASTA with motifs. 
    
    """
    import random
    
    all_records = []

    for uFS_id in joint_seqs:
        
        nbirths = int(weights[uFS_id])
        weight = min(nbirths, max_weight)

        # Pick nbirths random sequences to add to the motif file
        for i in range(weight):
            rec = random.choice(joint_seqs[uFS_id])
            all_records.append(rec)
            
    if len(all_records) > max_total_seqs:
        log(f"[motif] subsampling to {max_total_seqs:,} sequences")
        all_records = random.sample(all_records, max_total_seqs)
            
    return all_records


def shuffle_foreground_seqs(foreground_records):
    """
    Instead of testing against a background of random sequences, 
    test against a background of shuffled foreground sequences in order
    to preserve the length and base composition.

    Args:
        foreground_records (_type_): _description_
    """
    # Shuffle each foreground sequence to preserve length + composition
    fg_shuffled = []
    for fg_rec in foreground_records:
        seq_list = list(str(fg_rec.seq))
        random.shuffle(seq_list)
        fg_seq = ''.join(seq_list)
        fg_shuffled.append(SeqRecord(Seq(fg_seq), id=f"shuffled_{fg_rec.id}"))       
    return fg_shuffled

def build_background(genome_seq, lengths, seed):
    """
    Sample random genome windows for background model. 
    Use the median length of the flanking sequences.
    
    """
    random.seed(seed)
    n      = len(genome_seq)
    records = []
    n_samples = len(lengths)
    #attempts = 0
    #max_attempts = n_samples * 100

    for length in lengths:
        #while len(records) < n_samples and attempts < max_attempts:
        #    attempts += 1
        start = random.randint(0, n - length)
        seq = genome_seq[start:start + length]
        
        if len(seq) != length or "N" in seq:
            continue
        
        rec = SeqRecord(Seq(seq), id=f"bg_{len(records)}", name="", description=f'{start}-{start + length}')
        records.append(rec)

    log(f"[background] sampled {len(records):,}  (target {n_samples}, attempts x)")
    return records


# ── main ─────────────────────────────────────────────────────────────────────

log("=== motif_prep.py starting ===")

flanking_sequences = load_flanking_sequences(snakemake.input.anchors_tsv)
births_5 = load_births(snakemake.input.births_5)
births_3 = load_births(snakemake.input.births_3)
births = {**births_5, **births_3}
genome_seq = load_reference(snakemake.input.reference)

# Join 5' and 3' flanking sequences
joint_seqs = join_flanking_sequences(
    snakemake.input.reference_insertions, 
    flanking_sequences, 
    snakemake.input.anchor_map,
    min_length = min_length,
    max_length = max_length
)

# Write to fasta
weighted_records = write_motif_fasta(
    joint_seqs, births, 
    max_weight = max_weights, 
    max_total_seqs = max_total_seqs
    )

SeqIO.write(weighted_records, snakemake.output.fg, "fasta")
log(f"[write] {snakemake.output.fg}  ({len(weighted_records):,} seqs)")


# Get background sequences: same number and lengths as foreground
foreground_lengths = [len(rec.seq) for rec in weighted_records]

# Background (use median size of foreground sequences, randomly sampled from genome)
for i, bg_path in enumerate(snakemake.output.bg):
    bg_records = build_background(
        genome_seq, foreground_lengths,
        seed = snakemake.params.seed + i   # different seed per file
    )
    SeqIO.write(bg_records, bg_path, "fasta")
    log(f"[write] {bg_path}  ({len(bg_records):,} seqs)")


log("=== motif_prep.py done ===")
log_fh.close()