#!/bin/bash -l

#SBATCH --job-name=data_analysis_abcd
#SBATCH --output=data_analysis_abcd_%A_%a.out
#SBATCH --error=data_analysis_abcd_%A_%a.err

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16g
#SBATCH -t 24:00:00

#SBATCH --account=mfiecas
#SBATCH --mail-type=ALL
#SBATCH --mail-user=neher015@umn.edu

cd /home/mfiecas/neher015/abcd_multiview

module load R

Rscript 03_data_analysis_abcd.R