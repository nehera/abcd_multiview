#!/bin/bash -l
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16g
#SBATCH --tmp=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=neher015@umn.edu
module load R/4.0.4
R CMD BATCH mcmc-v5.R mcmc-v5.codeoutput.txt