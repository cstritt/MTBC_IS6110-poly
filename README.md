# Evolutionary dynamics of IS6110 in the *Mycobacterium tuberculosis* complex

This folder contains the material of the study "A niche constraints model for the evolutionary dynamics of the insertion sequence IS6110 in the *Mycobacterium tuberculosis* complex", including the Snakemake workflow to replicate the full analysis. 

The results are available as a bioRxiv preprint: https://doi.org/10.64898/2026.09.24.753865.

The tool used to infer IS polymorphisms in the study, detettore6110, is available here: https://github.com/cstritt/detettore6110. 

## Snakemake workflow

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