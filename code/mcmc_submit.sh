#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=16g
#SBATCH -t 4:00:00
#SBATCH --account=mfiecas
#SBATCH --mail-type=ALL
#SBATCH --mail-user=neher015@umn.edu

cd /home/mfiecas/neher015/abcd_multiview

module load R

Rscript code/mcmc_master.R $SLURM_ARRAY_TASK_ID