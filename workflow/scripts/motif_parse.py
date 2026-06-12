import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist

def parse_motif_pwm(xml_path):
    """Extract all motifs and their PWMs from a STREME XML file."""
    tree = ET.parse(xml_path)
    root = tree.getroot()
    motifs = []
    for motif in root.iter("motif"):
        pwm = []
        for pos in motif:
            if pos.tag == "pos":
                pwm.append([float(pos.get("A")), float(pos.get("C")),
                            float(pos.get("G")), float(pos.get("T"))])
        motifs.append({
            "id":        motif.get("id"),
            "consensus": motif.get("alt"),
            "width":     int(motif.get("width")),
            "pvalue":    float(motif.get("train_pvalue")),
            "pwm":       np.array(pwm)   # shape (width, 4)
        })
    return motifs

def pwm_distance(pwm_a, pwm_b):
    """
    Compute minimum Euclidean distance between two PWMs,
    allowing for offset alignment and reverse complement.
    Returns the minimum distance across all alignments.
    """
    def rc_pwm(pwm):
        return pwm[::-1, [3, 2, 1, 0]]   # reverse rows, swap A<->T C<->G

    def align_distance(a, b):
        # Slide shorter over longer, return minimum column-wise RMSE
        if len(a) > len(b):
            a, b = b, a
        w = len(a)
        dists = []
        for start in range(len(b) - w + 1):
            d = np.sqrt(np.mean((a - b[start:start+w])**2))
            dists.append(d)
        return min(dists)

    d_fwd = align_distance(pwm_a, pwm_b)
    d_rc  = align_distance(rc_pwm(pwm_a), pwm_b)
    return min(d_fwd, d_rc)

# ── Load all background XMLs ──────────────────────────────────────────────────
all_motifs = []
for bg_i, xml_path in enumerate(snakemake.input.xml_bg):
    for motif in parse_motif_pwm(xml_path):
        motif["background"] = bg_i
        all_motifs.append(motif)


# ── Cluster motifs by PWM similarity ─────────────────────────────────────────
# Compute pairwise distances
n = len(all_motifs)
dist_matrix = np.zeros((n, n))
for i in range(n):
    for j in range(i+1, n):
        d = pwm_distance(all_motifs[i]["pwm"], all_motifs[j]["pwm"])
        dist_matrix[i, j] = d
        dist_matrix[j, i] = d

# Hierarchical clustering
Z = linkage(pdist(dist_matrix), method="average")
labels = fcluster(Z, t=0.3, criterion="distance")   # threshold controls merging sensitivity

# ── Summarise clusters ────────────────────────────────────────────────────────
summary = []
for cluster_id in np.unique(labels):
    idx = np.where(labels == cluster_id)[0]
    cluster_motifs = [all_motifs[i] for i in idx]

    n_backgrounds = len(set(m["background"] for m in cluster_motifs))
    consensuses   = [m["consensus"] for m in cluster_motifs]
    pvalues       = [m["pvalue"]    for m in cluster_motifs]

    # Average PWM (align to shortest first)
    min_width = min(m["pwm"].shape[0] for m in cluster_motifs)
    pwm_stack = np.stack([m["pwm"][:min_width] for m in cluster_motifs])
    mean_pwm  = pwm_stack.mean(axis=0)

    summary.append({
        "cluster_id":     cluster_id,
        "n_backgrounds":  n_backgrounds,        # how many of 10 runs found this motif
        "n_motifs":       len(cluster_motifs),
        "consensuses":    ";".join(set(consensuses)),
        "median_pvalue":  np.median(pvalues),
        "mean_pwm":       mean_pwm
    })

summary_df = pd.DataFrame(summary).sort_values("n_backgrounds", ascending=False)

# ── Write output ──────────────────────────────────────────────────────────────
summary_df.drop(columns="mean_pwm").to_csv(snakemake.output.tsv, sep="\t", index=False)
