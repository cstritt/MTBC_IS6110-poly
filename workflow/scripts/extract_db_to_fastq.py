#!/usr/bin/env python3
"""
Extract sequences from SQLite database and write to gzipped FASTQ file.
This script extracts all sequences for a given strain from the IS6110 partial hits database.
"""

import sqlite3
import sys
import gzip
from pathlib import Path


def extract_strain_to_fastq(db_path, strain_name, output_path):
    """
    Extract all sequences for a given strain from the database and write to gzipped FASTQ.
    
    Args:
        db_path: Path to the SQLite database
        strain_name: Sample/strain name to extract
        output_path: Path for the output gzipped FASTQ file
    """
    # Connect to database with timeout for concurrent access
    conn = sqlite3.connect(db_path, timeout=60)
    conn.execute("PRAGMA journal_mode=WAL;")
    cursor = conn.cursor()
    
    # Query all sequences for this strain
    query = """
    SELECT sequence_key, fasta_sequence
    FROM sequences
    WHERE sample_name = ?
    ORDER BY id
    """
    
    cursor.execute(query, (strain_name,))
    sequences = cursor.fetchall()
    
    # Write to gzipped FASTQ
    # Use a fixed high-quality score (I = Phred 40) since we don't have original qualities
    with gzip.open(output_path, 'wt') as f:
        for seq_id, seq_data in sequences:
            # Write FASTQ format
            f.write(f"@{seq_id}\n")
            f.write(f"{seq_data}\n")
            f.write(f"+\n")
            # Use fixed quality score (I = Phred 40) for all bases
            f.write("I" * len(seq_data) + "\n")
    
    conn.close()
    
    return len(sequences)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("Usage: extract_db_to_fastq.py <database_path> <strain_name> <output_path>")
    
    db_path = sys.argv[1]
    strain_name = sys.argv[2]
    output_path = sys.argv[3]
    
    # Create output directory if it doesn't exist
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    
    count = extract_strain_to_fastq(db_path, strain_name, output_path)
    print(f"Extracted {count} sequences for strain {strain_name} to {output_path}")
