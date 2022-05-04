#!/bin/bash
#SBATCH --job-name=fl_MC_lj
#SBATCH --output=sbatchout/fl_MC_lj%a.out
#SBATCH --error=sbatchout/fl_MC_lj%a.err
#SBATCH --ntasks=1
#SBATCH --array=0-29
#SBATCH --mem=15G
#SBATCH --partition=hns,normal

LAMBDA_LIST=($(<../../lambda_leggauss.dat))
lambda=${LAMBDA_LIST[${SLURM_ARRAY_TASK_ID}]}
echo "running simulation at lambda" $lambda
../../out/build/hs $lambda $SLURM_ARRAY_TASK_ID phase1
