#!/bin/bash

#SBATCH --mail-type ALL 
#SBATCH --mail-user Senbai.Kang@chuv.ch

#SBATCH --output out

#SBATCH --partition cpu

#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=500gb

module load singularity

singularity exec --bind /work,/scratch /scratch/skang/containers/latest_bayesspace.sif Rscript --slave debug.R
