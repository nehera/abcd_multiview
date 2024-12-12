#!/bin/bash -l

#SBATCH --job-name=simulation_study
#SBATCH --output=tmp/simulation_study_%A_%a.out
#SBATCH --error=tmp/simulation_study_%A_%a.err

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH -t 01:30:00

#SBATCH --account=mfiecas
#SBATCH --mail-type=ALL
#SBATCH --mail-user=neher015@umn.edu

cd $HOME/abcd_multiview

module load R/4.3.0-openblas-rocky8

Rscript simulation_study.R $SLURM_ARRAY_TASK_ID
