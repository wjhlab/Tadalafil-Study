#!/bin/bash
#$ -N Tadalafil_
#$ -cwd
#$ -pe local 4
#$ -l h_vmem=64G,mem_free=64G,h_fsize=50G
#$ -m e
#$ -M zzhan220@jh.edu
#$ -e ./report/
#$ -o ./report/

module load singularity/3.2.1
module list

DATA_DIR="loom_files/"
LOOM_FILES=`find ${DATA_DIR} | grep "[[:digit:]]_counts.loom$"`
PREFIX="loom_files/tadalafil_"
SUFFIX="_counts.loom"

# HG38 refrence files downloaded from aertslab/github

mkdir "processed_data/05_pyscenic"

LOOM="loom_files/tadalafil_VEH_counts.loom";
PT_ID=${LOOM#"$PREFIX"}
PT_ID=${PT_ID%"$SUFFIX"}
  
RESULT_DIR="processed_data/05_pyscenic/J1568_${PT_ID}"
GRN="${RESULT_DIR}/adj_${PT_ID}.tsv"
  
if [ -f "${GRN}" ]
then 
  echo "${PT_ID} already has a grn result. Proceeding to next subject"
else  
  echo "${PT_ID} : Calculating gene regulatory network"
  mkdir "${RESULT_DIR}"
    
  singularity exec ~/aertslab-pyscenic-0.11.0.sif pyscenic grn \
	${LOOM} \
   	mm_mgi_tfs.txt \
    	-o $GRN \
    	--num_workers 4 \
        --seed 123

  echo "${PT_ID} Step one: adjacency matrix calculation complete"
fi

exit 0
