library(Seurat)
library(domino)

source("DominoScript/class_definitions.R")
source("DominoScript/import_fxns.R")

sessionInfo()

# list directories with pyscenic results
pyscenic_path <- "processed_data/05_pyscenic"
pyscenic_res <- list.dirs(path = pyscenic_path, recursive = FALSE)

result_dir <- "processed_data/06_domino"
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}

# download cellphoneDB reference data
cell_db_path <- "data/cell_db"
if(!dir.exists(cell_db_path)){
  dir.create(cell_db_path, recursive = TRUE)
}

# check if the cellphone db files have been downloaded
complexes_check <- file.exists(paste0(cell_db_path, "/complexes.csv"))
interactions_check <- file.exists(paste0(cell_db_path, "/interactions.csv"))
proteins_check <- file.exists(paste0(cell_db_path, "/proteins.csv"))
genes_check <- file.exists(paste0(cell_db_path, "/genes.csv"))

if(!complexes_check){
  system(
    'curl -O https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/sources/complex_curated.csv \
    mv complex_curated.csv data/cell_db/complexes.csv')}
if(!interactions_check){
  system(
    'curl -O https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/sources/interaction_curated.csv \
  mv interaction_curated.csv data/cell_db/interactions.csv')}
if(!proteins_check){
  system(
    'curl -O https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/sources/protein_curated.csv \
  mv protein_curated.csv data/cell_db/proteins.csv')}
if(!genes_check){
  system(
    'curl -O https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/gene_input.csv \
  mv gene_input.csv data/cell_db/genes.csv')}

# cellphoneDB files downloaded 1:00 pm 6/16/22
# cellphoneDB V2.0.0

res <- "processed_data/05_pyscenic/CLIST"
pt_id <- "CLIST"
print(res)
print(pt_id)

cat(paste0(pt_id), " : Loading Pyscenic Results")

# Patient's Seurat object
seurat <- readRDS("processed_data/Seurat_CLIST.rds")

# Pyscenic results data
auc <- read.table(paste0(res, "/auc_", pt_id, ".csv"),
                  header = TRUE, row.names = 1,
                  stringsAsFactors = FALSE, sep = ",")
regulons <- paste0(res, "/regulons_", pt_id, ".csv")

# seurat object information
counts <- seurat@assays$RNA@counts
z_scores <- as.matrix(seurat@assays$RNA@scale.data)
clusters <- as.factor(seurat$cell_type)

# Create domino object
domino <- create_domino(signaling_db = "data/cell_db",
                        features = t(auc),
                        counts = counts,
                        z_scores = z_scores,
                        clusters = clusters,
                        df = regulons,
                        gene_conv = c("HGNC", "MGI"),
                        gene_conv_host = "dec2021.archive.ensembl.org",
                        remove_rec_dropout = FALSE)

# save the unbuilt domino object for changing build parameters
saveRDS(domino, file = paste0(result_dir, "/", pt_id, "_domino_unbuilt.rds"))

cat(paste0(pt_id), " : Domino Complete")


# clear global environment to prevent the workspace from being saved
rm(list = ls())

