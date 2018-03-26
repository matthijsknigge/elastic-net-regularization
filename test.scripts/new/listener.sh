#!/bin/bash
#SBATCH --job-name=listener
#SBATCH --output=listener.out
#SBATCH --err=listener.err
#SBATCH --time=00:50:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --mem=25gb

module load Python/3.5.1-foss-2015b

python3 listener.py

