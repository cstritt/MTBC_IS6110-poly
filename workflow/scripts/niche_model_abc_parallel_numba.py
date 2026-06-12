#!/usr/bin/env python3
"""
Approximate Bayesian Computation for niche model with parallel execution.

Priors:
  - p_essential: Uniform[0.5, 0.95]  (most sites are essential)
  - p_tolerated: Uniform[0.05, 0.5]  (minority of non-essential are tolerated)
  - r_birth: Uniform[0.001, 2.0]     (per-copy transposition rate)
  - r_purge: Uniform[0.001, 10.0]    (purging rate on tolerated sites)
"""

import os
import sys
import random
import numpy as np
import pandas as pd
from multiprocessing import Pool, cpu_count
from ete3 import Tree

import niche_model_numba as niche_model


def run_one_simulation(args):
    """
    Worker function for parallel ABC.
    
    Args:
        args: tuple of (sim_id, tree, L, tip_names, lineage_idx, unique_lineages)
    
    Returns:
        (sim_id, params_list, stats_dict)
    """
    sim_id, tree_str, L, tip_names, lineage_idx, unique_lineages = args
    
    # Seed each worker independently based on PID and sim_id
    worker_seed = (os.getpid() * 1000003 + sim_id * 97) % (2**31)
    np.random.seed(worker_seed)
    random.seed(worker_seed)
    
    # Sample from priors
    p_essential = random.uniform(0.4, 0.95)
    p_tolerated = random.uniform(0.2, 0.8)
    r_birth = random.uniform(0.00001, 0.0002)
    r_purge = random.uniform(0.0001, 1)
    
    # Parse tree (reconstruct from Newick string)
    tree = Tree(tree_str, format=1)
    
    # Run simulation
    try:
        params, cn, stats = niche_model.run_simulation(
            tree=tree,
            L=L,
            p_essential=p_essential,
            p_tolerated=p_tolerated,
            r_birth=r_birth,
            r_purge=r_purge,
            tip_names=tip_names,
            lineage_idx=lineage_idx,
            unique_lineages=unique_lineages,
            seed=worker_seed
        )
        
        params_list = [p_essential, p_tolerated, r_birth, r_purge]
        
        return (sim_id, params_list, stats)
    
    except Exception as e:
        print(f"Simulation {sim_id} failed: {e}", file=sys.stderr)
        return (sim_id, [np.nan]*4, {})


def get_worker_count():
    """Get correct number of workers respecting SLURM allocation."""
    # Check if running under SLURM
    if 'SLURM_CPUS_PER_TASK' in os.environ:
        return int(os.environ['SLURM_CPUS_PER_TASK'])
    elif 'SLURM_NTASKS' in os.environ:
        return int(os.environ['SLURM_NTASKS'])
    else:
        # Fallback for local/non-SLURM execution
        return cpu_count()
    

def main():
    """
    Main ABC pipeline.
    
    Usage:
        python niche_model_abc_parallel_numba.py <tree.nwk> <nsim> <metadata.tsv> <outdir>
    """
    try:
        treepath = sys.argv[1]
        nsim = int(sys.argv[2])
        metadata_path = sys.argv[3]
        outdir = sys.argv[4]
    except IndexError:
        print("Usage: niche_model_abc_parallel_numba.py <tree.nwk> <nsim> <metadata.tsv> <outdir>",
              file=sys.stderr)
        sys.exit(1)
    
    # Ensure output directory exists
    os.makedirs(outdir, exist_ok=True)
    
    # Load tree
    tree = Tree(treepath, format=1)
    tree_str = tree.write(format=1)  # Convert to Newick for serialization
    tip_names = sorted(tree.get_leaf_names())
    
    # Load metadata to map tips to lineages
    metadata = pd.read_csv(metadata_path, sep='\t')
    
    # Create lineage index
    lineage_map = dict(zip(metadata['GNUMBER'], metadata['LINEAGE_x']))
    unique_lineages = sorted(set(lineage_map.values()))
    lineage_idx = np.array([unique_lineages.index(lineage_map.get(tip, 'unknown'))
                            for tip in tip_names], dtype=int)
    
    # Genome length (from MTBC0)
    L = 7178
    
    # Get number of workers
    n_workers = get_worker_count()    
    
    print(f"Tree: {len(tip_names)} tips")
    print(f"Genome length: {L}")
    print(f"Lineages: {unique_lineages}")
    print(f"Running {nsim} simulations with {n_workers} workers...\n")
    
    # Prepare task arguments
    tasks = [
        (i, tree_str, L, tip_names, lineage_idx, unique_lineages)
        for i in range(nsim)
    ]
    
    # Run in parallel
    params_list = []
    stats_list = []
    
    with Pool(processes=n_workers) as pool:
        for i, (sim_id, params, stats) in enumerate(pool.imap_unordered(run_one_simulation, tasks)):
            if i % 100 == 0:
                print(f"Completed {i}/{nsim} simulations")
            
            params_list.append(params)
            stats_list.append(stats)
    
    print(f"\nCompleted all {nsim} simulations\n")
    
    # Assemble dataframes
    params_df = pd.DataFrame(
        params_list,
        columns=['p_essential', 'p_tolerated', 'r_birth', 'r_purge']
    )
    
    # All stats dicts should have same keys; build uniform dataframe
    if stats_list and stats_list[0]:
        stats_df = pd.DataFrame(stats_list)
    else:
        stats_df = pd.DataFrame()
    
    # Write outputs
    params_path = os.path.join(outdir, 'abc_params.tsv')
    stats_path = os.path.join(outdir, 'abc_summaries.tsv')
    
    params_df.to_csv(params_path, sep='\t', index=False)
    stats_df.to_csv(stats_path, sep='\t', index=False)
    
    print(f"Saved {params_path}")
    print(f"Saved {stats_path}")
    
    # Summary
    print("\nParameter summary:")
    print(params_df.describe())
    
    print("\nStatistics summary:")
    print(stats_df.describe())


if __name__ == "__main__":
    main()
