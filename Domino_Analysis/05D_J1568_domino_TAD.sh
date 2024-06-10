#!/bin/bash
#$ -N Tadalafil_domino_
#$ -cwd
#$ -pe local 1
#$ -l h_vmem=36G,mem_free=36G,h_fsize=50G
#$ -m e
#$ -M zzhan220@jh.edu
#$ -e ./report/
#$ -o ./report/

module load conda_R/4.1.x
module list

R_SCRIPT="scripts/06_J1568_domino_TAD.R"

R CMD BATCH $R_SCRIPT

echo "Domino Results Calculation Complete"

exit 0
