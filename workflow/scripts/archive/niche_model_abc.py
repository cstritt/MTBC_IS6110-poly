#!/usr/bin/env python3

"""
Approximate Bayesian Computation

"""

import niche_model
import numpy as np
import os
import pandas as pd
import random
import sys

from ete3 import Tree
from scipy import stats

def get_summary_stats(copy_numbers: pd.DataFrame, metadata: pd.DataFrame):
    
    # Map lineage info from metadata using the first column as index
    copy_numbers = copy_numbers.copy()
    
    # Assume metadata has columns: ['SAMPLE_ID', 'LINEAGE_x']
    merged = copy_numbers.merge(metadata[['GNUMBER', 'LINEAGE_x']],
                                left_on=copy_numbers.columns[0],
                                right_on='GNUMBER',
                                how='left')
    
    # Summary statistics for the second column
    summary = stats.describe(merged.iloc[:, 1].values)
    
    # Group by lineage and compute mean, std for the second column
    cladevar = merged.groupby('LINEAGE_x')[merged.columns[1]].agg(['mean', 'std'])
    cladevar.columns = ['Mean', 'SD']

    # Compute variance and std of lineage means
    cladevar_var = np.var(cladevar['Mean'], ddof=1)
    cladevar_sd = np.std(cladevar['Mean'], ddof=1)

    # Assemble summary data in a dictionary
    result = {
        'mean': summary.mean,
        'std': np.sqrt(summary.variance),
        'min': summary.minmax[0],
        '25%': np.percentile(copy_numbers.iloc[:, 1], 25),
        '50%': np.percentile(copy_numbers.iloc[:, 1], 50),
        '75%': np.percentile(copy_numbers.iloc[:, 1], 75),
        'max': summary.minmax[1],
        'cladevar': cladevar_var,
        'cladesd': cladevar_sd
    }

    return pd.Series(result)


def main():
    
    try:
        treepath = sys.argv[1]
        nsim = int(sys.argv[2])
        metadata = pd.read_csv(sys.argv[3], sep='\t')
        outpath = sys.argv[4]
    except IndexError:
        sys.exit('Usage: niche_model_abc.py <tree> <nsim> <metadata> <outpath>')    
    
    # Fixed parameters
    initial_copies = 1
    max_copies = 40

    seed = 42
    random.seed(seed)
    np.random.seed(seed)

    # Read tree with ete3
    tree = Tree(treepath, format=1)  # format=1 tries to preserve branch lengths and names
    root = tree.get_tree_root()         

    # Define empty matrices
    params_mat = np.empty(shape=(nsim, 3), dtype=float)
    params_mat[:] = np.nan
    summaries_mat = np.empty(shape=(nsim, 9), dtype=float)
    summaries_mat[:] = np.nan

    # Run simulations
    for i in range(nsim):
        
        # Sample from prior distributions
        p_lethal = round(random.uniform(0.1, 0.99), 2)
        r_birth = int(random.uniform(1000, 10000))
        r_death = int(random.uniform(100, 5000))

        params_mat[i] = [p_lethal, r_birth, r_death]

        params = {
            "p_lethal": p_lethal,
            "r_birth": r_birth,
            "r_death": r_death,
            "max_copies": max_copies
        }
        
        results = {}
        niche_model.traverse_and_simulate_ete(root, initial_copies, params, results)
        copy_numbers = pd.DataFrame([{"tip": k, "copies": v} for k, v in sorted(results.items())])
        
        cn_stats = get_summary_stats(copy_numbers, metadata)
        summaries_mat[i] = cn_stats
        
    # Save tables to files
    pd.DataFrame(params_mat).to_csv(os.path.join(outpath, 'abc_params.tsv'), index=False, sep='\t') 
    pd.DataFrame(summaries_mat).to_csv(os.path.join(outpath, 'abc_summaries.tsv'), index=False, sep='\t')

if __name__ == "__main__":
    main()