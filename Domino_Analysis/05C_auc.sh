#!/bin/bash
#$ -N Tadalafil_auc
#$ -cwd
#$ -pe local 8
#$ -l h_vmem=36G,mem_free=36G,h_fsize=50G
#$ -m e
#$ -M zzhan220@jh.edu
#$ -e ./report/
#$ -o ./report/

module load singularity/3.2.1
module list

DATA_DIR="processed_data/05_pyscenic"
PT_DIRS=$(ls -d ${DATA_DIR}/*)
PREFIX="processed_data/05_pyscenic/"

for PT in $PT_DIRS; do
  PT_ID=${PT#"$PREFIX"}
  REGULON="${PT}/regulons_${PT_ID}.csv"
  LOOM="loom_files/${PT_ID}_counts.loom"
  AUC="${PT}/auc_${PT_ID}.csv"
  
  if [ -f "${AUC}" ]; then 
    echo "${PT_ID} already has a regulon result. Proceeding to next subject"
  else 
    echo "${PT_ID} : Calculating AUC matrix"
    
    singularity exec ~/aertslab-pyscenic-0.11.0.sif pyscenic aucell \
  	"${LOOM}" \
  	"${REGULON}" \
  	-o "${AUC}"
    
    echo "${PT_ID} Step 3: AUC matrix calculation complete"
  fi
done

exit 0
