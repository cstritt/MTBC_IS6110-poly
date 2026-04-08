#!/bin/bash

#SBATCH --job-name={jobname}
#SBATCH --partition={resources.slurm_partition}
#SBATCH --qos={resources.qos}
#SBATCH --cpus-per-task={resources.cpus_per_task}
#SBATCH --mem={resources.mem_mb}
#SBATCH --time={resources.runtime}  # Must be SLURM format: D-HH:MM:SS or HH:MM:SS
#SBATCH --output=log/%x.%j.out                 # where to store the output ( %j is the jobID )
#SBATCH --error=log/%x.%j.err                  # where to store error messages (%x is the jobname)

#Run .bashrc to initialize conda and julia
source $HOME/.bashrc

# Activate conda env
# conda activate (env-name)
{exec_job}