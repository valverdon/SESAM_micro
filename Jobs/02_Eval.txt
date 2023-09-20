#!/bin/bash

#SBATCH --job-name=Eval$1$2$3
#SBATCH --account=aguisan_sometalp
#SBATCH --time=10:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --mem=10G
#SBATCH --array=1-251
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --mail-user valentin.verdon@unil.ch
#SBATCH --mail-type BEGIN,END,FAIL



source /dcsrsoft/spack/bin/setup_dcsrsoft
module load gcc/9.3.0 r/4.0.5
module load geos/3.8.1 
module load netcdf-c/4.8.0
module load proj/5.2.0 
module load gdal/2.4.4-proj-5.2.0

Rscript code/02_Eval.R $SLURM_ARRAY_TASK_ID $1 $2 $3