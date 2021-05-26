#!/bin/bash

#SBATCH -p athey
#SBATCH -c 1
#SBATCH -G 1
#SBATCH --mem=16G

module load python/3.6.1
module load py-pytorch/1.4.0_py36
module load py-numpy/1.19.2_py36
srun python3 generation/estimate_gans.py
