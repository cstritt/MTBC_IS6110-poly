"""
scripts/motif_parse.py
======================
Snakemake script — called by rule motif_parse.

Parses STREME XML output for both the 5' and 3' anchor sides and writes a
single summary TSV with one row per discovered motif.

Output columns
--------------
side                  5prime / 3prime
motif_id              STREME motif identifier (e.g. 1-ACGT)
consensus             IUPAC consensus string (alt attribute in XML)
width                 motif width (bp)
log10_pvalue          log10(p-value) from discriminative test (negative = significant)
pvalue                raw p-value
train_pos_count       foreground sequences containing the motif
train_neg_count       background sequences containing the motif
dtc                   discriminative total count
score_threshold       position-weight-matrix score threshold used
"""

import xml.etree.ElementTree as ET
import pandas as pd

log_path = str(snakemake.log[0])
log_fh   = open(log_path, "w")

def log(msg):
    log_fh.write(msg + "\n")
    log_fh.flush()


def parse_streme_xml(xml_path):
    """
    Parse a STREME streme.xml file.
    Returns a list of dicts, one per motif.
    """
    log(f"[parse] {xml_path}  ")
    tree = ET.parse(xml_path)
    root = tree.getroot()

    motifs = []
    for motif in root.iter("motif"):
        a = motif.attrib

        # p-value: STREME reports train_log10pvalue (negative, more negative = better)
        raw_log10p = a.get("train_log10pvalue", "NA")
        try:
            log10p = float(raw_log10p)
        except ValueError:
            log10p = float("nan")

        raw_p = a.get("train_pvalue", "NA")
        try:
            pval = float(raw_p)
        except ValueError:
            pval = float("nan")

        motifs.append({
            "motif_id":         a.get("id",                  ""),
            "consensus":        a.get("alt",                 ""),   # IUPAC-like
            "width":            a.get("width",               ""),
            "log10_pvalue":     log10p,
            "pvalue":           pval,
            "train_pos_count":  a.get("train_pos_count",     ""),
            "train_neg_count":  a.get("train_neg_count",     ""),
            "dtc":              a.get("train_dtc",           ""),
            "score_threshold":  a.get("score_threshold",     ""),
        })

    log(f"[parse] found {len(motifs)} motifs")
    return motifs


log("=== motif_parse.py starting ===")

motifs = parse_streme_xml(snakemake.input.xml)
motifs_bg = parse_streme_xml(snakemake.input.xml_bg)

# Add column specificying bg or no bg model, join data frames
for m in motifs:
    m["test"] = 'shuffle'
for m in motifs_bg:
    m["test"] = 'background'
    
motifs.extend(motifs_bg)

if motifs:
    df = pd.DataFrame(motifs)
    # Sort: most significant first (most negative log10p)
    df = df.sort_values("log10_pvalue")
    df.to_csv(snakemake.output.tsv, sep="\t", index=False)
    log(f"[write] {snakemake.output.tsv}  ({len(df)} motifs total)")
else:
    # Write empty file so Snakemake output is satisfied
    pd.DataFrame(columns=[
        "side", "motif_id", "consensus", "width",
        "log10_pvalue", "pvalue",
        "train_pos_count", "train_neg_count", "dtc", "score_threshold"
    ]).to_csv(snakemake.output.tsv, sep="\t", index=False)
    log("[write] no motifs found – wrote empty TSV")

log("=== motif_parse.py done ===")
log_fh.close()
