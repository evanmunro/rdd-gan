#!/bin/bash

#SBATCH -p athey
#SBATCH -c 1
#SBATCH -G 1
#SBATCH --mem=16G

module load python/3.9.0
source rddgan/bin/activate
srun python generation/estimate_gans.py
