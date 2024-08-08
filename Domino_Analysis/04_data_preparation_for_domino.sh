#!/bin/bash
#$ -N prep_pyscenic_domino
#$ -cwd
#$ -pe local 2
#$ -l h_vmem=30G,mem_free=30G,h_fsize=50G
#$ -m e
#$ -M jmitch81@jhmi.edu
#$ -e ./report/
#$ -o ./report/

module load conda_R/4.1.x
module list

R_SCRIPT="scripts/04_data_preparation_for_domino.R"

pwd

echo "Begining pyscenic preprocessing"

R CMD BATCH $R_SCRIPT

echo "Preprocessing complete"

exit 0
