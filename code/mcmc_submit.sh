#!/bin/bash -l

#SBATCH --job-name=mcmc_submit
#SBATCH --array=1-5
#SBATCH --output=mcmc_submit_%A_%a.out
#SBATCH --error=mcmc_submit_%A_%a.err

#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --mem=16g
#SBATCH -t 48:00:00

#SBATCH --account=mfiecas
#SBATCH --mail-type=ALL
#SBATCH --mail-user=neher015@umn.edu

cd /home/mfiecas/neher015/abcd_multiview

module load R

Rscript code/mcmc_master.R $SLURM_ARRAY_TASK_ID