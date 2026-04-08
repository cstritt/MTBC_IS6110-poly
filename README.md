# Transposition and decay if IS in the MTBC

This repository contains a complete Snakemake workflow to replicate the analysis of IS6110 in x MTBC strains, 
as described here: . 

Getting input data and running the workflow are split

# Get data
data/get_data.ipynb

# Set up snakemake



# Run main workflow
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

# Provide alternative config file
snakemake \
  --profile cluster \
  --configfile config/config_georgia.yaml \
  --dry-run \
  main



```

# Run detettore benchmarking
```{bash}

snakemake --profile cluster -n benchmarking
snakemake --profile cluster benchmarking


```

# Run detettore only with alternative config
```{bash}

snakemake --profile cluster -n detettore_only --configfile config/config_tanzania.yml
snakemake --profile cluster detettore_only


```

# Run workflow on personal computer
```{bash}
snakemake benchmarking -n
snakemake benchmarking --cores 8 --rerun-incomplete --use-conda
```