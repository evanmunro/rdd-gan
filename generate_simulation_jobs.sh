#!/bin/bash

#SBATCH -p athey
#SBATCH -c 1
#SBATCH -t 08:00:00
#SBATCH --mem=16G

module load R/4.0.2
Rscript simulations/slurm_simulation.R mc_curved
Rscript simulations/slurm_simulation.R curved
