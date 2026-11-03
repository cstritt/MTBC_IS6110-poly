#!/usr/bin/env python3
"""
parse_gff_regions.py

Extracts a region table from a PGAP-annotated GFF3 file.
Optionally annotates each region with overlapping insertion sequences
from an ISEScan GFF file and computes GC content from a reference FASTA.

Output columns
--------------
    region_name   gene name / locus tag, or "geneA;geneB" for intergenic gaps
    region_type   feature type (gene, rRNA, tRNA, pseudogene…) or "intergenic"
    chromosome    sequence/contig name (GFF column 1)
    start         1-based start coordinate
    end           1-based end coordinate
    length        end - start + 1
    gc_content    GC fraction (0–1), ignoring non-ACGT bases
    strand        +, -, or . (intergenic / unknown)
    gene_id       ID attribute from column 9 (PGAP locus tag)
    gene_name     Name/gene attribute, if present
    product       product description, if present
    pseudo        "true" if annotated as pseudogene, else "false"
    note          Note attribute, if present
    is_overlap    "true" / "false" — whether the region overlaps any IS element
    is_ids        semicolon-separated IS IDs that overlap
    is_family     semicolon-separated IS families
    is_cluster    semicolon-separated IS clusters
    is_overlap_bp total base pairs of overlap with IS elements
"""

import argparse
import csv
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_attributes(attr_str: str) -> dict:
    """Parse a GFF3 column-9 attribute string into a key->value dict."""
    attrs = {}
    for field in attr_str.rstrip(";").split(";"):
        field = field.strip()
        if "=" in field:
            key, _, value = field.partition("=")
            attrs[key.strip()] = value.strip()
    return attrs


def safe_get(attrs: dict, *keys, default=".") -> str:
    """Return the value of the first matching key, or default."""
    for k in keys:
        if k in attrs:
            return attrs[k]
    return default


def intergenic_name(locus_a: str, locus_b: str) -> str:
    return f"{locus_a};{locus_b}"


def overlap_bp(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """Return overlapping bp between two 1-based inclusive intervals."""
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def read_fasta(path: str) -> dict[str, str]:
    """Read FASTA into dict: contig -> uppercase sequence."""
    seqs = {}
    name = None
    chunks = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

        if name is not None:
            seqs[name] = "".join(chunks).upper()

    return seqs


def calc_gc(seq: str) -> float | str:
    """Calculate GC fraction, ignoring non-ACGT bases."""
    seq = seq.upper()

    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")

    denom = gc + at
    if denom == 0:
        return "."

    return gc / denom


def annotate_gc(row: dict, genome: dict[str, str]) -> dict:
    """Annotate row with gc_content."""
    chrom = row["chromosome"]

    if chrom not in genome:
        row["gc_content"] = "."
        return row

    seq = genome[chrom]

    start = row["start"] - 1
    end = row["end"]

    if start < 0 or end > len(seq):
        row["gc_content"] = "."
        return row

    subseq = seq[start:end]
    gc = calc_gc(subseq)

    if gc == ".":
        row["gc_content"] = "."
    else:
        row["gc_content"] = f"{gc:.4f}"

    return row


# ---------------------------------------------------------------------------
# ISEScan parsing
# ---------------------------------------------------------------------------

def parse_isescan(isescan_path: str) -> list[dict]:
    elements = []

    with open(isescan_path) as fh:
        for line in fh:
            line = line.rstrip("\n")

            if line.startswith("#") or not line.strip():
                continue

            cols = line.split("\t")
            if len(cols) < 9:
                continue

            seqname, source, ftype, start, end, score, strand, phase, attr_str = cols[:9]

            if ftype != "insertion_sequence":
                continue

            attrs = parse_attributes(attr_str)

            elements.append({
                "seqname": seqname,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "is_id": attrs.get("ID", "."),
                "family": attrs.get("family", "."),
                "cluster": attrs.get("cluster", "."),
            })

    return elements


def build_is_index(is_elements: list[dict]) -> dict:
    index = defaultdict(list)

    for el in is_elements:
        index[el["seqname"]].append(el)

    return dict(index)


def annotate_is_overlap(row: dict, is_index: dict) -> dict:
    chrom = row["chromosome"]
    r_start = row["start"]
    r_end = row["end"]

    overlapping = []

    for el in is_index.get(chrom, []):
        bp = overlap_bp(r_start, r_end, el["start"], el["end"])
        if bp > 0:
            overlapping.append((el, bp))

    if overlapping:
        row["is_overlap"] = "true"
        row["is_ids"] = ";".join(el["is_id"] for el, _ in overlapping)
        row["is_family"] = ";".join(el["family"] for el, _ in overlapping)
        row["is_cluster"] = ";".join(el["cluster"] for el, _ in overlapping)
        row["is_overlap_bp"] = sum(bp for _, bp in overlapping)
    else:
        row["is_overlap"] = "false"
        row["is_ids"] = "."
        row["is_family"] = "."
        row["is_cluster"] = "."
        row["is_overlap_bp"] = 0

    return row


# ---------------------------------------------------------------------------
# GFF parsing
# ---------------------------------------------------------------------------

def parse_gff(gff_path: str, feature_types: set) -> tuple[dict, dict]:
    features_by_seq = defaultdict(list)
    sequence_regions = {}

    with open(gff_path) as fh:
        for line in fh:
            line = line.rstrip("\n")

            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    sequence_regions[parts[1]] = int(parts[3])
                continue

            if line.startswith("#") or not line.strip():
                continue

            cols = line.split("\t")
            if len(cols) < 9:
                continue

            seqname, source, ftype, start, end, score, strand, phase, attr_str = cols[:9]

            if feature_types and ftype not in feature_types:
                continue

            attrs = parse_attributes(attr_str)

            feat = {
                "chromosome": seqname,
                "region_type": ftype,
                "start": int(start),
                "end": int(end),
                "strand": strand if strand in ("+", "-") else ".",
                "gene_id": safe_get(attrs, "ID", "locus_tag"),
                "gene_name": safe_get(attrs, "Name", "gene"),
                "product": safe_get(attrs, "product"),
                "locus_tag": safe_get(attrs, "locus_tag"),
                "note": safe_get(attrs, "Note", "note"),
                "pseudo": "true" if ("pseudo" in attrs or ftype == "pseudogene") else "false",
            }

            features_by_seq[seqname].append(feat)

    for seq in features_by_seq:
        features_by_seq[seq].sort(key=lambda x: (x["start"], x["end"]))

    return dict(features_by_seq), sequence_regions


# ---------------------------------------------------------------------------
# Table construction
# ---------------------------------------------------------------------------

def finalize_row(row: dict, is_index: dict | None, genome: dict[str, str] | None) -> dict:
    if is_index is not None:
        row = annotate_is_overlap(row, is_index)

    if genome is not None:
        row = annotate_gc(row, genome)
    else:
        row["gc_content"] = "."

    return row


def build_table(
    features_by_seq: dict,
    sequence_regions: dict,
    genome_length: int | None,
    include_intergenic: bool,
    is_index: dict | None,
    genome: dict[str, str] | None,
) -> list[dict]:

    rows = []

    for seqname, feats in sorted(features_by_seq.items()):
        chrom_len = sequence_regions.get(seqname, genome_length)

        prev_end = 0
        prev_locus = "start"

        for feat in feats:

            if include_intergenic and feat["start"] > prev_end + 1:
                ig_start = prev_end + 1
                ig_end = feat["start"] - 1

                curr_locus = (
                    feat["locus_tag"]
                    if feat["locus_tag"] != "."
                    else feat["gene_id"]
                )

                row = {
                    "region_name": intergenic_name(prev_locus, curr_locus),
                    "region_type": "intergenic",
                    "chromosome": seqname,
                    "start": ig_start,
                    "end": ig_end,
                    "length": ig_end - ig_start + 1,
                    "strand": ".",
                    "gene_id": ".",
                    "gene_name": ".",
                    "product": ".",
                    "pseudo": ".",
                    "note": ".",
                }

                rows.append(finalize_row(row, is_index, genome))

            row = {
                "region_name": feat["locus_tag"],
                "region_type": feat["region_type"],
                "chromosome": feat["chromosome"],
                "start": feat["start"],
                "end": feat["end"],
                "length": feat["end"] - feat["start"] + 1,
                "strand": feat["strand"],
                "gene_id": feat["gene_id"],
                "gene_name": feat["gene_name"],
                "product": feat["product"],
                "pseudo": feat["pseudo"],
                "note": feat["note"],
            }

            rows.append(finalize_row(row, is_index, genome))

            locus = feat["locus_tag"] if feat["locus_tag"] != "." else feat["gene_id"]

            if feat["end"] > prev_end:
                prev_end = feat["end"]
                prev_locus = locus

        if include_intergenic and chrom_len and prev_end < chrom_len:
            ig_start = prev_end + 1
            ig_end = chrom_len

            row = {
                "region_name": intergenic_name(prev_locus, "end"),
                "region_type": "intergenic",
                "chromosome": seqname,
                "start": ig_start,
                "end": ig_end,
                "length": ig_end - ig_start + 1,
                "strand": ".",
                "gene_id": ".",
                "gene_name": ".",
                "product": ".",
                "pseudo": ".",
                "note": ".",
            }

            rows.append(finalize_row(row, is_index, genome))

    return rows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

BASE_COLUMNS = [
    "region_name",
    "region_type",
    "chromosome",
    "start",
    "end",
    "length",
    "gc_content",
    "strand",
    "gene_id",
    "gene_name",
    "product",
    "pseudo",
    "note",
]

IS_COLUMNS = [
    "is_overlap",
    "is_ids",
    "is_family",
    "is_cluster",
    "is_overlap_bp",
]

DEFAULT_FEATURE_TYPES = {
    "gene",
    "pseudogene",
    "rRNA",
    "tRNA",
    "ncRNA",
    "tmRNA",
    "repeat_region",
    "regulatory",
}


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--gff", required=True)
    parser.add_argument("--out", default="regions.tsv")
    parser.add_argument("--isescan", default=None)
    parser.add_argument("--fasta", default=None, help="Reference genome FASTA for GC annotation")
    parser.add_argument("--feature_types", nargs="*", default=None)
    parser.add_argument("--intergenic", dest="intergenic", action="store_true", default=True)
    parser.add_argument("--no_intergenic", dest="intergenic", action="store_false")
    parser.add_argument("--genome_length", type=int, default=None)

    args = parser.parse_args()

    feature_types = set(args.feature_types) if args.feature_types else DEFAULT_FEATURE_TYPES

    print(f"[INFO] Input GFF       : {args.gff}", file=sys.stderr)
    print(f"[INFO] ISEScan GFF     : {args.isescan or 'not provided'}", file=sys.stderr)
    print(f"[INFO] Reference FASTA : {args.fasta or 'not provided'}", file=sys.stderr)

    is_index = None
    if args.isescan:
        is_elements = parse_isescan(args.isescan)
        is_index = build_is_index(is_elements)
        print(f"[INFO] IS elements     : {len(is_elements)} loaded", file=sys.stderr)

    genome = None
    if args.fasta:
        genome = read_fasta(args.fasta)
        print(f"[INFO] FASTA contigs    : {len(genome)} loaded", file=sys.stderr)

    features_by_seq, sequence_regions = parse_gff(args.gff, feature_types)

    rows = build_table(
        features_by_seq=features_by_seq,
        sequence_regions=sequence_regions,
        genome_length=args.genome_length,
        include_intergenic=args.intergenic,
        is_index=is_index,
        genome=genome,
    )

    columns = BASE_COLUMNS + (IS_COLUMNS if is_index is not None else [])

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"[INFO] Wrote {len(rows)} rows → {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()

