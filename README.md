# Evolutionary dynamics of IS6110 in the Mycobacterium tuberculosis complex

This repository contains the workflow used in the study "Evolutionary dynamics of the insertion sequence IS6110 in the *Mycobacterium tuberculosis* complex". It is organized as a Snakemake workflow for quantifying IS6110 insertion polymorphisms and benchmarking detection performance. 

The preprint is on bioRxiv: https://doi.org/10.64898/2026.09.24.753865.


## Repository structure

- `workflow/` — Snakemake pipeline and cluster configuration
- `data/` — reference genomes, annotations, variant files
- `metadata/` — sample metadata and curated lookup tables used by the analysis
- `figures/` — generated manuscript and supplementary figures
- `manuscript/` — manuscripts, drafts, notes, and figure planning documents

The study considers 10,000 strains such that full replication would require a lot of computational time. To provide a short cut, the following genomic data is provided through Zenodo (work in progress):

  - a VCF file with SNPs for all strains, which was used to infer the phylogenetic tree
  - a SQLite database containing only the short reads of all strains that partially map to IS6110. These are the reads used to detect IS6110 polymorphisms. 


## Data sources and inputs

The workflow is configured in `workflow/config/config.yaml` and supports three input modes:

1. `scicore` (default): uses existing CRAM files on the sciCORE HPC system
2. `db`: extracts IS6110 partial hits from a SQLite database ()

The relevant configuration includes paths to:

- reference genomes
- sample metadata
- annotated genomic regions
- benchmarking datasets
- IS6110 reference sequence

## Reproducing the analysis

### 1) Create the Snakemake environment

```bash
conda create -n snakemake -c conda-forge \
  snakemake snakemake-executor-plugin-slurm snakemake-storage-plugin-fs
```

### 2) Run a dry run on a cluster using SLURM

```bash
conda activate snakemake
screen -R snakemake

snakemake \
  --profile cluster \
  --dry-run
```

### 3) Run the full workflow on a cluster

```bash
snakemake --profile cluster
```

### 4) Run the workflow locally

```bash
snakemake --cores 8 \
  --rerun-incomplete \
  --use-conda
```

## Main workflow outputs

The main analysis entry point includes outputs for:

- benchmarking summaries
- copy-number analyses
- birth-death model summaries
- niche-space analyses
- motif-enrichment results
- figures and summary reports


## Related resources

- `detettore6110`: https://github.com/cstritt/detettore6110
- `large_variable_alignment`: https://github.com/cstritt/large_variable_alignment
