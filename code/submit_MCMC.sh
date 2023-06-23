#!/bin/bash -l
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH -p amdlarge
#SBATCH --ntasks=1
#SBATCH --mem=100g
#SBATCH --tmp=400g

#SBATCH --account=carr0603
#SBATCH --mail-type=ALL
#SBATCH --mail-user=neher015@umn.edu

cd /home/carr0603/neher015/Documents/GitHub/sim.khv
module load R/4.0.4

id=$SLURM_ARRAY_TASK_ID
R CMD BATCH "--args i=$id" code/05_run.sim.viz.R outs/reports/sim.viz$id.Rout

# this script is meant to be submitted in a job array (https://www.msi.umn.edu/support/faq/how-do-i-use-job-array)
# e.g. sbatch --array=1-1000 jobtemplate_2.sh
