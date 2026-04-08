#!/usr/bin/env python3
import argparse
import random
from collections import Counter
import numpy as np
from scipy.stats import chisquare
import pandas as pd
import matplotlib.pyplot as plt

# Optional logo
try:
    import logomaker
    HAVE_LOGOMAKER = True
except ImportError:
    HAVE_LOGOMAKER = False


def read_fasta(path):
    seqs = []
    with open(path) as f:
        seq = ""
        for line in f:
            line = line.strip()
            if not line: 
                continue
            if line.startswith(">"):
                if seq:
                    seqs.append(seq.upper())
                    seq = ""
            else:
                seq += line
        if seq:
            seqs.append(seq.upper())
    return seqs


def shuffle_seq(seq):
    s = list(seq)
    random.shuffle(s)
    return "".join(s)


def trim_center(seq, w):
    """Trim to w bp centered on the middle."""
    mid = len(seq) // 2
    half = w // 2
    start = mid - half
    end = start + w
    return seq[start:end]


def build_pfm(seqs):
    """Return a dictionary {A:[...],C:[...],G:[...],T:[...]}."""
    L = len(seqs[0])
    pfm = {b: [0]*L for b in "ACGT"}
    for s in seqs:
        for i, ch in enumerate(s):
            if ch in pfm:
                pfm[ch][i] += 1
    return pfm


def freq_matrix(pfm):
    """Convert counts to frequencies."""
    L = len(pfm["A"])
    freq = {b: [] for b in "ACGT"}
    for i in range(L):
        total = sum(pfm[b][i] for b in "ACGT")
        for b in "ACGT":
            freq[b].append(pfm[b][i] / total if total > 0 else 0)
    return freq


def chi_square_positions(pos_seqs, bg_seqs):
    L = len(pos_seqs[0])
    pvals = []
    for i in range(L):
        pos_counts = Counter(s[i] for s in pos_seqs)
        bg_counts  = Counter(s[i] for s in bg_seqs)

        pos_vec = np.array([pos_counts[b] for b in "ACGT"])
        bg_vec  = np.array([bg_counts[b] for b in "ACGT"])

        # Expected = background frequencies scaled to positives count
        bg_freq = bg_vec / bg_vec.sum()
        expected = bg_freq * pos_vec.sum()

        chi2, p = chisquare(pos_vec, expected)
        pvals.append(p)
    return pvals


def make_logo(freq, outfile):
    if not HAVE_LOGOMAKER:
        print("logomaker not installed; skipping logo.")
        return
    df = pd.DataFrame(freq)
    fig, ax = plt.subplots(figsize=(12,2))
    logo = logomaker.Logo(df, ax=ax)
    ax.set_ylabel("frequency")
    ax.set_xlabel("Position")
    plt.tight_layout()
    plt.savefig(outfile, dpi=200)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Test sequence preference near IS insertion sites.")
    parser.add_argument("--f5", required=True, help="5' flank FASTA")
    parser.add_argument("--f3", required=True, help="3' flank FASTA")
    parser.add_argument("-w", "--window", type=int, default=10, help="window size (bp)")
    args = parser.parse_args()

    # Load sequences
    f5 = read_fasta(args.f5)
    f3 = read_fasta(args.f3)
    seqs = f5 + f3

    # Trim to window
    seqs = [trim_center(s, args.window) for s in seqs]

    # Background: shuffle each sequence independently
    bg = [shuffle_seq(s) for s in seqs]

    # Build PFM and freq matrices
    pfm_pos = build_pfm(seqs)
    freq_pos = freq_matrix(pfm_pos)

    pfm_bg = build_pfm(bg)
    freq_bg = freq_matrix(pfm_bg)

    # Position-wise chi-square test
    pvals = chi_square_positions(seqs, bg)

    # Save statistics
    df = pd.DataFrame({
        "position": range(args.window),
        "pval": pvals,
        "A_freq_pos": freq_pos["A"],
        "C_freq_pos": freq_pos["C"],
        "G_freq_pos": freq_pos["G"],
        "T_freq_pos": freq_pos["T"],
        "A_freq_bg": freq_bg["A"],
        "C_freq_bg": freq_bg["C"],
        "G_freq_bg": freq_bg["G"],
        "T_freq_bg": freq_bg["T"],
    })
    df.to_csv("position_stats.tsv", sep="\t", index=False)
    print("Saved: position_stats.tsv")

    # Logo of the real flanks
    make_logo(freq_pos, "sequence_logo.png")
    print("Saved: sequence_logo.png")


if __name__ == "__main__":
    main()
