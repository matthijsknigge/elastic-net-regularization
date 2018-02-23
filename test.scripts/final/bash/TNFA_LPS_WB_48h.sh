
#!/bin/bash
#SBATCH --job-name=TNFA_LPS_WB_48h
#SBATCH --output=log/TNFA_LPS_WB_48h.out
#SBATCH --err=log/TNFA_LPS_WB_48h.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --mem=20gb

module load R

Rscript cytokines.R -c TNFA_LPS_WB_48h
			   