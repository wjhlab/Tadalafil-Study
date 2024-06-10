# 04_data_preparation_for_domino
# Jacob Mitchell
# 12/06/2022

library(Seurat)
library(SeuratDisk)
library(hdf5r)
library(loomR)
library(dplyr)
library(ggplot2)

sessionInfo()

# directory for saved figures
figurepath <- "E:/WJHLab/domino/figure"
if(!dir.exists(figurepath)){
  dir.create(figurepath, recursive = TRUE)
}
# directory for files 
resultpath <- "E:/WJHLab/domino/result"
if(!dir.exists(resultpath)){
  dir.create(resultpath, recursive = TRUE)
}
# directory for cell count reports
countpath <- paste0(resultpath, "/counts_report")
if(!dir.exists(countpath)){
  dir.create(countpath, recursive = TRUE)
}

# read in seurat object
seurat <- readRDS("E:/WJHLab/domino/filtered_sobj_withSampleid.RDS")
head(rownames(seurat@assays$RNA@counts))
# report cell counts before processing
report_counts <- function(seurat, var.1, var.2 = NULL, filename){
  cfile_prefix <- "cell_count_"
  if(is.null(var.2)){
    writeLines(
      capture.output(table(seurat@meta.data[[var.1]])),
      con = paste0(countpath, "/", cfile_prefix, filename, ".txt")
    )
  } else {
    writeLines(
      capture.output(table(seurat@meta.data[[var.1]], 
                           seurat@meta.data[[var.2]])),
      con = paste0(countpath, "/", cfile_prefix, filename, ".txt")
    )
  }
}

# cell counts per arm
report_counts(seurat = seurat,
              var.1 = "group_id",
              filename = "group_id")
report_counts(seurat = seurat,
              var.1 = "Arm",
              var.2 = "cell_type",
              filename = "all_Arm-cell_type")
report_counts(seurat = seurat,
              var.1 = "Arm",
              var.2 = "PubID",
              filename = "all_Arm-PubID")

# cell counts per OS group
report_counts(seurat = seurat,
              var.1 = "OS",
              filename = "all_OS")
report_counts(seurat = seurat,
              var.1 = "OS",
              var.2 = "cell_type",
              filename = "all_OS-cell_type")
report_counts(seurat = seurat,
              var.1 = "OS",
              var.2 = "PubID",
              filename = "all_OS-PubID")

# Scaling of expression values will occur per patient

# save loom file of counts data
# all cells
if(!file.exists(paste0(resultpath, "/all_counts.loom"))){
  loom_all <- create(filename = paste0(resultpath, "/all_counts.loom"),
                     data = seurat@assays$RNA@counts)
  loom_all$close_all()
}

# rds is not saved because no changes have been made since the 03 result

# subset seurat object by each patient
pt_ser_list <- list()
pt_ids <- unique(seurat@meta.data$group_id)


pt_ser_list[["CLIST"]] <- seurat[,seurat$group_id == "CLIST"]

# check cell counts for correct subsetting
report_counts(seurat = pt_ser_list[["CLIST"]], 
              var.1 = "cell_type",
              filename = paste0("CLIST", "_CellType"))

# Scale expression values
pt_ser_list[["CLIST"]] <- ScaleData(pt_ser_list[["CLIST"]],
                               features = rownames(pt_ser_list[["CLIST"]]))

# Save loom file of counts data
if(!file.exists(paste0(resultpath, "/tadalafil_", "CLIST", "_counts.loom"))){
  loom <- create(filename = paste0(resultpath, "/tadalafil_", "CLIST", "_counts.loom"),
                 data = pt_ser_list[["CLIST"]]@assays$RNA@counts)
  loom$close_all()
}

# Save seurat object as rds
saveRDS(object = pt_ser_list[["CLIST"]],
        file = paste0(resultpath, "/J1568_Seurat_", "CLIST", ".rds"))


ser_veh = subset(seurat,subset = group_id == "VEH")
a = list("Gene" = row.names(ser_veh@assays$RNA@counts))
lfile$add.col.attribute(attribute = a, overwrite = TRUE)
lfile[['row_attrs']]
maat[1:5,1:5]

lfile$close_all()

lfile <- connect(filename = paste0("E:/WJHLab/domino/",
                                   "result/",
                                   "tadalafil_VEH_counts.loom"),
                 mode = "r+")
# matrix needs to be transposed as loom files default to cell x feature
maat <- (lfile[["matrix"]][,])
maat[1:5,1:5]

class(colnames(ser_veh@assays$RNA@counts))


