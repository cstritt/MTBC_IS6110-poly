#%%
"""
Aim: get a fasta with all features and intergenic regions in MTBC0, 
plus a description in the header with the ortholog/homolog in H37Rv.

Ortholog assignment strategy (three-tier):
  1. Reciprocal Best Hit (RBH): MTBC0->H37Rv best hit must also hit back the 
     same MTBC0 gene. High confidence, used for most unique genes.
  2. Synteny tiebreaker: for repetitive families (PPE, PE-PGRS, IS elements)
     where multiple MTBC0 genes map equally well to the same H37Rv gene,
     the candidate whose normalised chromosomal rank is closest to the MTBC0
     gene's rank wins.
  3. Unassigned: genes with no BLAST hit above thresholds get 'NA'.

The 'method' column in the output records how each assignment was made.

GFF structure notes (from GFFExaminer):
  MTBC0 (PGAP):
    - 'gene'       (4012): parent for all protein-coding genes; gene_biotype qualifier
                           distinguishes protein_coding / tRNA / rRNA / tmRNA etc.
    - 'pseudogene' ( 146): standalone parent rows for pseudogenes
    - tRNA/rRNA/CDS/exon: child rows — skipped here to avoid double-counting

  H37Rv (RefSeq):
    - 'gene'       (3978): parent for protein-coding genes + some RNA genes
    - 'pseudogene' (  30): standalone (no 'gene' parent)
    - 'tRNA'       (  45): standalone (no 'gene' parent in RefSeq format)
    - 'rRNA'       (   3): standalone
    - 'ncRNA'      (  20): standalone
    => All five types must be included to capture the full H37Rv feature set.
"""

import os
import sys
import argparse
import subprocess
from collections import defaultdict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from BCBio import GFF


def parse_args():
    parser = argparse.ArgumentParser(description='Create MTBC0 to H37Rv gene map')
    parser.add_argument('--mtbc0_annot', required=True, help='MTBC0 GFF annotation file')
    parser.add_argument('--mtbc0_genome', required=True, help='MTBC0 genome FASTA')
    parser.add_argument('--h37rv_annot', required=True, help='H37Rv GFF annotation file')
    parser.add_argument('--h37rv_genome', required=True, help='H37Rv genome FASTA')
    parser.add_argument('--output', required=True, help='Output genemap TSV file')
    return parser.parse_args()

args = parse_args()

mtbc0_annot  = args.mtbc0_annot
mtbc0_genome = args.mtbc0_genome
h37rv_annot  = args.h37rv_annot
h37rv_genome = args.h37rv_genome
outpath     = args.output

# ── feature type filters (one per GFF, based on GFFExaminer output) ──────────

# MTBC0 PGAP: 'gene' captures all biotypes (protein_coding, tRNA, rRNA, tmRNA,
# ncRNA ...) via the gene_biotype qualifier. 'pseudogene' is a separate parent.
MTBC0_GENE_TYPES = ['gene', 'pseudogene']

# H37Rv RefSeq: 'gene' captures protein-coding + some RNA genes. tRNA, rRNA,
# ncRNA, and pseudogene are standalone entries without a 'gene' parent.
H37RV_GENE_TYPES = ['gene', 'pseudogene', 'tRNA', 'rRNA', 'ncRNA']

# ── thresholds ────────────────────────────────────────────────────────────────

MIN_COVERAGE  = 0.50   # minimum query coverage fraction to keep a hit
MIN_PIDENT    = 70.0   # minimum % identity to keep a hit
BITSCORE_TOL  = 0.02   # hits within 2 % of top bitscore are treated as tied


# ── BLAST wrapper ─────────────────────────────────────────────────────────────

class Blast:
    def __init__(self):
        # query_id -> list of hit dicts, sorted by bitscore descending
        self.results = {}

    def seqs_to_files(self, query, target):
        SeqIO.write(query,  'query.fa',  'fasta')
        SeqIO.write(target, 'target.fa', 'fasta')

    def blastn(self, max_target_seqs=10):
        """
        Run blastn and store ALL hits per query (up to max_target_seqs).
        Hits are filtered by MIN_COVERAGE and MIN_PIDENT and sorted by
        bitscore descending.
        """
        cmd = [
            'blastn',
            '-query',           'query.fa',
            '-subject',         'target.fa',
            '-outfmt',          '6 qseqid qstart qend qlen sseqid sstart send pident length bitscore',
            '-max_target_seqs', str(max_target_seqs),
        ]

        proc   = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        output = proc.stdout.read()

        raw = defaultdict(list)
        for line in output.splitlines():
            line   = line.decode('utf-8')
            fields = line.split('\t')
            qseqid, qstart, qend, qlen, sseqid, sstart, send, pident, length, bitscore = fields
            qlen      = int(qlen)
            length    = int(length)
            pident    = float(pident)
            bitscore  = float(bitscore)
            qcov      = length / qlen

            if qcov < MIN_COVERAGE or pident < MIN_PIDENT:
                continue

            raw[qseqid].append({
                'sseqid':   sseqid,
                'pident':   pident,
                'length':   length,
                'qlen':     qlen,
                'bitscore': bitscore,
                'qcov':     qcov,
            })

        for q in raw:
            raw[q].sort(key=lambda h: -h['bitscore'])

        self.results = dict(raw)

    def clean_up(self):
        for f in ('query.fa', 'target.fa'):
            if os.path.exists(f):
                os.remove(f)


# ── helpers ───────────────────────────────────────────────────────────────────

def _get_feature_id(feature) -> str:
    """
    Return the most stable identifier for a feature.
    Priority: locus_tag > ID attribute.
    locus_tag is preferred because it is stable across annotation versions
    and is the key used in the genemap output.
    """
    lt = feature.qualifiers.get('locus_tag', [])
    if lt:
        return lt[0]
    return feature.qualifiers.get('ID', ['unknown'])[0]


def _get_gene_biotype(feature) -> str:
    """
    Return a human-readable biotype string.
    For PGAP 'gene' features this comes from gene_biotype qualifier.
    For H37Rv standalone tRNA/rRNA/ncRNA rows the feature.type itself is used.
    """
    bt = feature.qualifiers.get('gene_biotype', [])
    if bt:
        return bt[0]
    return feature.type   # e.g. 'tRNA', 'rRNA', 'pseudogene'


# ── feature extraction ────────────────────────────────────────────────────────

def extract_features(fasta_file, gff_file, limit_info):
    """
    Extract gene sequences and intergenic regions from genome + GFF.

    Parses only the top-level feature types specified in limit_info, so
    child rows (CDS, exon, tRNA child entries ...) are never seen and there
    is no double-counting.

    Returns:
      feature_d    : dict feature_id -> GFF feature object
      feature_seqs : list of SeqRecord (genes + intergenic, genes first)
      rank_d       : dict feature_id -> chromosomal rank (0-based, by position)
    """
    genome   = SeqIO.read(fasta_file, "fasta")
    inhandle = open(gff_file)

    feature_d    = {}
    feature_seqs = []
    feature_ids  = []
    feature_locs = []

    for rec in GFF.parse(inhandle, limit_info=limit_info):
        for feature in rec.features:

            start = int(feature.location.start)
            end   = int(feature.location.end)
            seq   = genome.seq[start:end]

            feature_id = _get_feature_id(feature)
            biotype    = _get_gene_biotype(feature)
            gene_name  = feature.qualifiers.get('gene', [''])[0]

            # Skip zero-length features (edge case in some annotations)
            if len(seq) == 0:
                continue

            # Avoid duplicates — defensive against limit_info edge cases
            if feature_id in feature_d:
                continue

            feature_rec = SeqRecord(
                seq, id=feature_id, name='',
                description=f'{gene_name}|{biotype}'
            )

            feature_d[feature_id] = feature
            feature_seqs.append(feature_rec)
            feature_ids.append(feature_id)
            feature_locs.append((start, end))

    inhandle.close()

    # Sort by genomic start position
    order               = sorted(range(len(feature_locs)), key=lambda i: feature_locs[i][0])
    feature_ids_sorted  = [feature_ids[i]  for i in order]
    feature_locs_sorted = [feature_locs[i] for i in order]
    feature_seqs_sorted = [feature_seqs[i] for i in order]

    # Intergenic regions between consecutive features
    intergenic_seqs = []
    for i in range(len(feature_locs_sorted) - 1):
        ig_start = feature_locs_sorted[i][1]
        ig_end   = feature_locs_sorted[i + 1][0]
        if ig_end > ig_start:
            seq   = genome.seq[ig_start:ig_end]
            ig_id = f'{feature_ids_sorted[i]}_{feature_ids_sorted[i + 1]}'
            intergenic_seqs.append(SeqRecord(seq, id=ig_id, name='', description=''))

    # Terminal regions
    if feature_locs_sorted:
        before = genome.seq[:feature_locs_sorted[0][0]]
        after  = genome.seq[feature_locs_sorted[-1][1]:]
        if len(before) > 0:
            intergenic_seqs.append(
                SeqRecord(before, id=f'0_{feature_ids_sorted[0]}', name='', description=''))
        if len(after) > 0:
            intergenic_seqs.append(
                SeqRecord(after, id=f'{feature_ids_sorted[-1]}_{len(genome)}',
                          name='', description=''))

    all_seqs = feature_seqs_sorted + intergenic_seqs
    rank_d   = {gid: rank for rank, gid in enumerate(feature_ids_sorted)}

    print(f'  {len(feature_d)} features extracted ({len(intergenic_seqs)} intergenic regions)')
    return feature_d, all_seqs, rank_d


# ── ortholog assignment ───────────────────────────────────────────────────────

def reciprocal_best_hits(fwd_blast: Blast, rev_blast: Blast) -> dict:
    """
    Returns dict: mtbc0_id -> h37rv_id for confirmed RBH pairs.
    """
    rbh = {}
    for mtbc0_id, hits in fwd_blast.results.items():
        best_h37rv = hits[0]['sseqid']
        rev_hits   = rev_blast.results.get(best_h37rv, [])
        if rev_hits and rev_hits[0]['sseqid'] == mtbc0_id:
            rbh[mtbc0_id] = best_h37rv
    return rbh


def synteny_tiebreak(mtbc0_id: str,
                     candidates: list,
                     mtbc0_ranks: dict,
                     h37rv_ranks: dict) -> str:
    """
    Among equally-scoring H37Rv candidates, pick the one whose normalised
    chromosomal rank is closest to the MTBC0 gene's rank.
    """
    n_mtbc0    = max(mtbc0_ranks.values()) + 1
    n_h37rv    = max(h37rv_ranks.values())  + 1
    mtbc0_norm = mtbc0_ranks.get(mtbc0_id, 0) / n_mtbc0

    best, best_dist = None, float('inf')
    for cand in candidates:
        h37rv_norm = h37rv_ranks.get(cand, 0) / n_h37rv
        dist = abs(mtbc0_norm - h37rv_norm)
        if dist < best_dist:
            best_dist = dist
            best = cand
    return best


def assign_orthologs(fwd_blast: Blast,
                     rev_blast: Blast,
                     feature_d: dict,
                     mtbc0_ranks: dict,
                     h37rv_ranks: dict) -> dict:
    """
    Three-tier assignment for every feature in feature_d:
      'RBH'       - reciprocal best hit (highest confidence)
      'synteny'   - synteny tiebreak among tied hits (repetitive families)
      'best_hit'  - unique best hit but not reciprocal (lower confidence)
      'no_hit'    - no BLAST hit above thresholds
    """
    rbh = reciprocal_best_hits(fwd_blast, rev_blast)

    assignments = {}
    for mtbc0_id in feature_d:

        if mtbc0_id in rbh:
            assignments[mtbc0_id] = (rbh[mtbc0_id], 'RBH')
            continue

        hits = fwd_blast.results.get(mtbc0_id, [])
        if not hits:
            assignments[mtbc0_id] = (None, 'no_hit')
            continue

        top_score  = hits[0]['bitscore']
        candidates = [h['sseqid'] for h in hits
                      if h['bitscore'] >= top_score * (1 - BITSCORE_TOL)]

        if len(candidates) == 1:
            assignments[mtbc0_id] = (candidates[0], 'best_hit')
        else:
            winner = synteny_tiebreak(mtbc0_id, candidates, mtbc0_ranks, h37rv_ranks)
            assignments[mtbc0_id] = (winner, 'synteny')

    return assignments


# ── MAIN ─────────────────────────────────────────────────────────────────────

#%% Extract features
print("Extracting MTBC0 features...")
feature_d, feature_seqs, mtbc0_ranks = extract_features(
    mtbc0_genome, mtbc0_annot,
    limit_info=dict(gff_type=MTBC0_GENE_TYPES)
)

print("Extracting H37Rv features...")
feature_d_h37rv, feature_seqs_h37rv, h37rv_ranks = extract_features(
    h37rv_genome, h37rv_annot,
    limit_info=dict(gff_type=H37RV_GENE_TYPES)
)

#%% Run BLAST in both directions (genes only, not intergenic)
gene_seqs_mtbc0 = [s for s in feature_seqs      if s.id in feature_d]
gene_seqs_h37rv = [s for s in feature_seqs_h37rv if s.id in feature_d_h37rv]

print(f"Running BLAST: MTBC0 -> H37Rv ({len(gene_seqs_mtbc0)} queries) ...")
mtbc0_vs_h37rv = Blast()
mtbc0_vs_h37rv.seqs_to_files(gene_seqs_mtbc0, gene_seqs_h37rv)
mtbc0_vs_h37rv.blastn(max_target_seqs=10)
mtbc0_vs_h37rv.clean_up()

print(f"Running BLAST: H37Rv -> MTBC0 ({len(gene_seqs_h37rv)} queries) ...")
h37rv_vs_mtbc0 = Blast()
h37rv_vs_mtbc0.seqs_to_files(gene_seqs_h37rv, gene_seqs_mtbc0)
h37rv_vs_mtbc0.blastn(max_target_seqs=10)
h37rv_vs_mtbc0.clean_up()

#%% Assign orthologs
print("Assigning orthologs...")
assignments = assign_orthologs(
    mtbc0_vs_h37rv, h37rv_vs_mtbc0,
    feature_d, mtbc0_ranks, h37rv_ranks,
)

from collections import Counter
method_counts = Counter(m for _, m in assignments.values())
print("Assignment summary:")
for method, count in sorted(method_counts.items()):
    print(f"  {method:15s}: {count}")

#%% Build rv_map and write output
rv_map     = {}
anc_only_l = []

for query, (h37rv_id, method) in assignments.items():

    feat    = feature_d[query]
    start   = int(feat.location.start)
    end     = int(feat.location.end)
    biotype = _get_gene_biotype(feat)

    if h37rv_id is None or h37rv_id not in feature_d_h37rv:
        anc_only_l.append(query)
        rv_map[query] = (start, end, biotype, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', method)
        continue

    rv_feat      = feature_d_h37rv[h37rv_id]
    rv_locus_tag = rv_feat.qualifiers.get('locus_tag', ['NA'])[0]
    rv_gene      = rv_feat.qualifiers.get('gene', [''])[0]
    rv_start     = int(rv_feat.location.start)
    rv_end       = int(rv_feat.location.end)
    rv_strand    = '+' if rv_feat.location.strand == 1 else '-'

    best_hit = mtbc0_vs_h37rv.results.get(query, [{}])[0]
    coverage = best_hit.get('qcov',   'NA')
    identity = best_hit.get('pident', 'NA')

    rv_map[query] = (start, end, biotype, rv_locus_tag, rv_start, rv_end,
                     rv_strand, rv_gene, coverage, identity, method)

print(f'Features with no H37Rv assignment: {len(anc_only_l)}')

#%% Write genemap TSV
header = ['mtbc0_gene_id', 'start', 'end', 'biotype',
          'h37rv_id', 'h37rv_start', 'h37rv_end', 'h37rv_strand',
          'gene_id', 'coverage', 'identity', 'method']

with open(outpath, 'w') as f:
    f.write('\t'.join(header) + '\n')
    for query, vals in rv_map.items():
        start, end, biotype, rv_id, rv_start, rv_end, rv_strand, gene_id, cov, ident, method = vals
        # locus_tag is already clean (no prefix stripping needed)
        f.write(f'{query}\t{start}\t{end}\t{biotype}\t{rv_id}\t{rv_start}\t{rv_end}\t'
                f'{rv_strand}\t{gene_id}\t{cov}\t{ident}\t{method}\n')

print(f"Written to {outpath}")