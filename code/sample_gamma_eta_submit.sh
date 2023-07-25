#!/bin/bash -l

#SBATCH --job-name=sample_gamma_eta
#SBATCH --array=1-3
#SBATCH --output=sample_gamma_eta_%A_%a.out
#SBATCH --error=sample_gamma_eta_%A_%a.err

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16g
#SBATCH -t 48:00:00

#SBATCH --account=mfiecas
#SBATCH --mail-type=ALL
#SBATCH --mail-user=neher015@umn.edu

cd /home/mfiecas/neher015/abcd_multiview

module load R

Rscript code/sample_gamma_eta.R $SLURM_ARRAY_TASK_ID