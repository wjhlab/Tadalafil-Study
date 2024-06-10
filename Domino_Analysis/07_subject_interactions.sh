#!/bin/bash
#$ -N Tadalafil_interactions_
#$ -cwd
#$ -pe local 8
#$ -l h_vmem=36G,mem_free=36G,h_fsize=50G
#$ -m e
#$ -M zzhan220@jh.edu
#$ -e ./report/
#$ -o ./report/

module load conda_R/4.1.x
module list

R_SCRIPT="scripts/07_subject_interactions.R"

echo "Compiling interactions from each subject's domino result"

R CMD BATCH $R_SCRIPT

echo "Compilation Complete"

exit 0
