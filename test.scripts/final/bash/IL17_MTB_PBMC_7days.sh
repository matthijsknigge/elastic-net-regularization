
#!/bin/bash
#SBATCH --job-name=IL17_MTB_PBMC_7days
#SBATCH --output=log/IL17_MTB_PBMC_7days.out
#SBATCH --err=log/IL17_MTB_PBMC_7days.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --mem=20gb

module load R

Rscript cytokines.R -c IL17_MTB_PBMC_7days
			   