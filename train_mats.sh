#!/bin/bash

#SBATCH -p athey
#SBATCH -c 1
#SBATCH -G 1
#SBATCH -t 08:00:00
#SBATCH --mem=64G

module load python/3.6.1
module load cuda/11.2.0
module load py-pytorch/1.4.0_py36
module load py-numpy/1.19.2_py36
module load py-pandas/1.0.3_py36
python3 generation/estimate_mats.py
