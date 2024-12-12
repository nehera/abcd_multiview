#!/bin/bash -l

#SBATCH --job-name=data_analysis
#SBATCH --output=tmp/data_analysis_%A_%a.out
#SBATCH --error=tmp/data_analysis_%A_%a.err

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH -t 96:00:00

#SBATCH --account=mfiecas
#SBATCH --mail-type=ALL
#SBATCH --mail-user=neher015@umn.edu

cd $HOME/abcd_multiview

module load R/4.3.0-openblas-rocky8

Rscript data_analysis.R $SLURM_ARRAY_TASK_ID