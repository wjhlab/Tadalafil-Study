#!/bin/bash
#$ -N Tadalafil_regulon_
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
  GRN="${PT}/adj_${PT_ID}.tsv"
  LOOM="loom_files/${PT_ID}_counts.loom"
  REGULON="${PT}/regulons_${PT_ID}.csv"
  
  if [ -f "${REGULON}" ]; then 
    echo "${PT_ID} already has a regulon result. Proceeding to next subject"
  else 
    echo "${PT_ID} : Calculating regulons"
    
singularity exec ~/aertslab-pyscenic-0.11.0.sif pyscenic ctx \
"${GRN}" \
mm10_motif_ref/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
mm10_motif_ref/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
--annotations_fname mm10_motif_ref/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl \
--expression_mtx_fname "${LOOM}" \
--mode "dask_multiprocessing" \
--output "${REGULON}" \
--num_workers 8
    
    echo "${PT_ID} Step 2: regulon calculation complete"
  fi
done

exit 0

  
