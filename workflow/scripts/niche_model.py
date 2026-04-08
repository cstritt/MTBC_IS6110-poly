#!/usr/bin/env python3
"""
simulate_is_with_niches_and_plot.py

Usage example:
    pip install ete3 matplotlib pandas numpy pillow
    python simulate_is_with_niches_and_plot.py --tree mytree.nwk \
      --p_lethal 0.25 --r_transpose 0.5 --r_delete 0.1 \
      --initial_copies 1 --max_niches 40 --output sim_results.csv \
      --fig fig_tree_hist.png --seed 42

Notes:
 - This script requires ete3 to render the tree. If you have rendering issues,
   ensure Pillow is installed: pip install pillow
 - The tree is rooted at midpoint if not already rooted.
"""

import argparse
import random
import numpy as np
import pandas as pd
import sys
from ete3 import Tree


def simulate_branch(start_copies, branch_length, r_transpose, r_delete, p_lethal, max_copies):
    """
    Simulates a single branch of a tree with IS elements.

    Args:
        start_copies (int): Number of copies of the IS element at the start of the branch.
        branch_length (float): Branch length in units of time.
        r_transpose (float): Rate of transposition events per copy per unit of time.
        r_delete (float): Rate of deletion events per copy per unit of time.
        p_lethal (float): Probability that a transposition event is lethal (i.e., it does not result in a successful insertion).
        max_copies (int): Maximum number of copies of the IS element allowed on the branch.

    Returns:
        int: The number of copies of the IS element at the end of the branch.
    """
    
    copies = int(start_copies)
    t = 0.0

    # Nothing can happen if branch length <= 0 or copies <= 0
    if branch_length <= 0.0 or copies <= 0:
        return copies

    while t < branch_length and copies > 0:
        per_copy_rate = r_transpose + r_delete
        total_rate = copies * per_copy_rate
        if total_rate <= 0.0:
            break

        wait = np.random.exponential(1.0 / total_rate)
        if t + wait > branch_length:
            break
        t += wait

        # Choose event type
        if random.random() < (r_transpose / per_copy_rate):
            # Transposition attempt
            if random.random() < p_lethal:
                # aborted insertion (no change)
                pass
            else:
                # Successful insertion: only increase if below max_niches
                if copies < max_copies:
                    copies += 1
                else:
                    # Option X: attempt occurs but has no effect
                    pass
        else:
            # Deletion event
            if copies > 0:
                copies -= 1
                # If copies becomes 0, births are impossible thereafter (absorbing 0)
    return copies


def simulate_branch_spatial(niches, branch_length, r_transpose, r_delete, p_lethal):
    """
    Spatially explicit IS model with new rule:
        - Transposition picks ANY genomic position uniformly.
        - Insertion only occurs if niches[pos] == 0 (tolerant + empty).
    
    niches:
        -1 = incompatible site (cannot host an insertion)
         0 = empty, insertion-compatible
         1 = occupied
    
    Returns:
        updated niches (copy)
    """

    niches = niches.copy()
    L = len(niches)

    t = 0.0

    while t < branch_length:
        # Where are the copies?
        occupied_indices = np.where(niches == 1)[0]
        copies = len(occupied_indices)

        # Absorbing state
        if copies == 0:
            return niches

        per_copy_rate = r_transpose + r_delete
        total_rate = copies * per_copy_rate
        if total_rate <= 0.0:
            return niches

        # Gillespie waiting time
        wait = np.random.exponential(1.0 / total_rate)
        if t + wait > branch_length:
            break
        t += wait

        # Choose event: transposition vs deletion
        if random.random() < (r_transpose / per_copy_rate):
            # --- Transposition attempt ---
            if random.random() < p_lethal:
                # aborted insertion
                continue

            # Choose ANY genomic position uniformly
            pos = random.randrange(L)

            # Check if insertion allowed
            if niches[pos] == 0:    # empty and tolerable
                niches[pos] = 1
            else:
                # incompatible or already occupied -> do nothing
                pass

        else:
            # --- Deletion event ---
            if copies > 0:
                pos = random.choice(occupied_indices)
                niches[pos] = 0

    return niches


def create_niche_space(L, p_lethal, initial_copies = 1):

    # Randomly distribut lethal niches
    """
    Create a vector of genomic niches compatible with IS insertions or not. 

    Parameters
    ----------
    L : int
        length of the niche space
    p_lethal : float
        probability of a site being lethal
    initial_copies : int, optional
        number of initial copies to seed in the niche space

    Returns
    -------
    niches : numpy array
        an array of length L with -1 (lethal), 0 (compatible) and 1 (occupied) values
    """
    niches = np.where(np.random.rand(L) < p_lethal, -1, 0) 

    # Seed initial copie(s)
    compat = np.where(niches == 0)[0]
    initial_sites = np.random.choice(compat, size=initial_copies, replace=False)
    niches[initial_sites] = 1

    return niches


def traverse_and_simulate_ete(node, inherited_copies, params, results):
    """
    Recursively simulate using an ete3 node.
    node.dist is the branch length from parent to this node.
    """
    branch_length = node.dist if node.dist is not None else 0.0
    
    end_copies = simulate_branch(
        inherited_copies,
        branch_length,
        params["r_birth"],
        params["r_death"],
        params["p_lethal"],
        params["max_copies"]
    )

    if node.is_leaf():
        name = node.name if node.name else "unnamed_tip"
        results[name] = end_copies
    else:
        for child in node.children:
            traverse_and_simulate_ete(child, end_copies, params, results)


def traverse_and_simulate_ete_spatial(node, inherited_niches, params, results):
    
    branch_length = node.dist if node.dist is not None else 0.0

    end_niches = simulate_branch_spatial(
        inherited_niches,
        branch_length,
        params["r_birth"],
        params["r_death"],
        params["p_lethal"]
    )

    if node.is_leaf():
        results[node.name or "unnamed_tip"] = end_niches
    else:
        for child in node.children:
            traverse_and_simulate_ete(child, end_niches, params, results)


def main():
    p = argparse.ArgumentParser(description="Simulate IS with niche cap and plot TREE | HISTOGRAM")
    p.add_argument("--tree", required=True, help="Input Newick tree file")
    p.add_argument("--p_lethal", type=float, default=0.95, help="Prob insertion hits lethal site (aborted)")
    p.add_argument("--r_birth", type=float, default=4000, help="Transposition rate per copy per unit branch length")
    p.add_argument("--r_death", type=float, default=0, help="Deletion rate per copy per unit branch length")
    p.add_argument("--initial_copies", type=int, default=1, help="Initial copies at root")
    p.add_argument("--max_copies", type=int, default=40, help="Maximum niches (cap) per lineage")
    p.add_argument("--spatial", action="store_true", default=False, help="Spatial simulation")
    p.add_argument("--niches", type=int, default=10000, help="Total number of niches, suitable or not")
    p.add_argument("--seed", type=int, default=42, help="Random seed")
    args = p.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)

    # Read tree with ete3
    tree = Tree(args.tree, format=1)  # format=1 tries to preserve branch lengths and names
    sys.stderr.write(tree.get_ascii(show_internal=True))

    params = {
        "p_lethal": args.p_lethal,
        "r_birth": args.r_birth,
        "r_death": args.r_death,
        "max_copies": args.max_copies
    }

    results = {}

    root = tree.get_tree_root()
    if args.spatial:
        niches = create_niche_space(args.niches, args.p_lethal, args.initial_copies)
        traverse_and_simulate_ete_spatial(root, niches, params, results)
    else:    
        traverse_and_simulate_ete(root, args.initial_copies, params, results)

    # Save results to CSV
    df = pd.DataFrame([{"tip": k, "copies": v} for k, v in sorted(results.items())])
    print(df.to_csv(sep="\t", index=False))


if __name__ == "__main__":
    main()