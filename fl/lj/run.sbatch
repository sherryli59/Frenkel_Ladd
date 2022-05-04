#!/bin/bash
#SBATCH --job-name=fl_MC_lj
#SBATCH --output=sbatchout/fl_MC_lj%a.out
#SBATCH --error=sbatchout/fl_MC_lj%a.err
#SBATCH --ntasks=1
#SBATCH --array=0-29
#SBATCH --mem=15G
#SBATCH --partition=hns,normal

LAMBDA_LIST=($(<../../lambda_leggauss.dat))
GAMMA_LIST=($(<../../gamma_leggauss.dat))
lambda=${LAMBDA_LIST[${SLURM_ARRAY_TASK_ID}]}
gamma=${GAMMA_LIST[${SLURM_ARRAY_TASK_ID}]}
echo "running phase1 simulation at lambda" $lambda
../../out/build/lj_phase1 $lambda $SLURM_ARRAY_TASK_ID phase1
echo "running phase2 simulation at gamma" $gamma
../../out/build/lj_phase2 $gamma $SLURM_ARRAY_TASK_ID phase2
