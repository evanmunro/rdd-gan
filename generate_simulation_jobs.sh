#!/bin/bash

#SBATCH -p athey
#SBATCH -c 1
#SBATCH -t 08:00:00
#SBATCH --mem=16G

module load R/4.0.2
Rscript simulation/slurm_simulation.R lee
Rscript simulation/slurm_simulation.R meyersson
Rscript simulation/slurm_simulation.R senate
Rscript simulation/slurm_simulation.R brazil
