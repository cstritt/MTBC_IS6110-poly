# Evolutionary dynamics of IS6110 in the *Mycobacterium tuberculosis* complex

This folder contains the material of the study "A niche constraints model for the evolutionary dynamics of the insertion sequence IS6110 in the *Mycobacterium tuberculosis* complex", including the Snakemake workflow to replicate the full analysis. 

The results are available as a bioRxiv preprint: https://doi.org/10.64898/2026.09.24.753865.

The tool used to infer IS polymorphisms in the study, detettore6110, is available here: https://github.com/cstritt/detettore6110. 

## Snakemake workflow

### Data source options

The workflow supports three data source options configured in `workflow/config/config.yaml`:

1. **scicore** (default): Uses existing CRAM files from sciCORE HPC storage
2. **ncbi**: Downloads FASTQ files from NCBI SRA using accession numbers from metadata
3. **db**: Extracts IS6110 partial hits from the SQLite database at `data/is6110/is6110_db/is6110.partial_hits.db`

To use the database source, set `data_source: db` in the config file.

### Set up snakemake

### Run main workflow
Install snakemake
```{bash}
conda create -n snakemake -c conda-forge snakemake snakemake-executor-plugin-slurm snakemake-storage-plugin-fs

```

Open screen, activate env, and run pipeline
```{bash}
screen -R snakemake 
conda activate snakemake

# Dry run
snakemake \
  --profile cluster \
  --dry-run \
  main

# Real run
snakemake --profile cluster main -n

```

### Run workflow on personal computer
```{bash}
snakemake benchmarking -n
snakemake benchmarking --cores 8 --rerun-incomplete --use-conda
```