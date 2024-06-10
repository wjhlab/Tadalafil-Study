library(domino)
library(Seurat)
library(ggplot2)
library(stringr)
library(scales)
library(grid)
library(circlize)
library(ComplexHeatmap)
library(dplyr)

setwd("E:/WJHLab/revised_domino/ligandList and signailing")

#Heatmap of TAD vs VEH for all cell type
########################################
sobj <- readRDS("E:/WJHLab/singlecell/filtered_sobj_withSampleid.RDS")
signaling1v2 <- readRDS("signaling1v2.RDS")
signaling1v2 <- signaling1v2[,order(colnames(signaling1v2))]
signaling4v5 <- readRDS("signaling4v5.RDS")
signaling4v5 <- signaling4v5[,order(colnames(signaling4v5))]
signalingTT <- readRDS("signalingTT.RDS")
signalingTT <- signalingTT[,order(colnames(signalingTT))]

sobj$cell_treat <- paste(sobj$cell_type, sobj$group_id, sep = "_")
cta <- colnames(signalingTT)
celltype <- gsub("_.*", "", cta)
treat_anno <- gsub(".*_", "", cta)
pal <- hue_pal()(length(unique(celltype)))
palette_key <- setNames(pal, unique(celltype))

pal_arm <- hue_pal()(length(unique(treat_anno)))
pal_arm_key <- setNames(pal_arm, unique(treat_anno))
column_ha <- HeatmapAnnotation(cell_type = celltype,
                               col = list(cell_type = palette_key))
arm_ha <- HeatmapAnnotation(Treatment = treat_anno,
                            col = list(treatment = pal_arm_key))

scaled1v2 <- as.matrix(t(scale(t(signalingTT))))
pdf("tadmlisttad_heatmap.pdf")
Heatmap(scaled1v2,
        name = "mean scaled\nexpression",
        col = colorRamp2(c(min(scaled1v2), 0, max(scaled1v2)),
                         c("blue", "#F9F9F9FF", "red")),
        top_annotation = arm_ha,
        bottom_annotation = column_ha,
        cluster_rows = TRUE, cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        column_names_gp = grid::gpar(fontsize = 4),
        row_names_gp = grid::gpar(fontsize = 4))
dev.off()


#Heatmap of MLIST vs MLIST.TAD for all cell type
################################################
signaling4v5 <- readRDS("signaling4v5.RDS")
signaling4v5 <- signaling4v5[,order(colnames(signaling4v5))]
signaling4v5 <- signaling4v5[rownames(signaling4v5) %in% c("Cxcl9", "Cxcl10", "Il12a", "Il12b", "Jag1", "Jag2"),]
signaling4v5 <- signaling4v5[, grepl("^(DC|Gran|Mac|Endo)", colnames(signaling4v5))]

sobj$cell_treat <- paste(sobj$cell_type, sobj$group_id, sep = "_")
cta <- colnames(signaling4v5)
celltype <- gsub("_.*", "", cta)
treat_anno <- gsub(".*_", "", cta)
pal <- hue_pal()(length(unique(celltype)))
palette_key <- setNames(pal, unique(celltype))

pal_arm <- hue_pal()(length(unique(treat_anno)))
pal_arm_key <- setNames(pal_arm, unique(treat_anno))
column_ha <- HeatmapAnnotation(cell_type = celltype,
                               col = list(cell_type = palette_key))
arm_ha <- HeatmapAnnotation(Treatment = treat_anno,
                            col = list(treatment = pal_arm_key))

scaled4v5 <- as.matrix(t(scale(t(signaling4v5))))
scaled4v5 <- na.omit(scaled4v5)
pdf("signaling4v5_heatmap_4type.pdf",height = 3, width = 5)
Heatmap(scaled4v5,
        name = "mean scaled\nexpression",
        col = colorRamp2(c(min(scaled4v5), 0, max(scaled4v5)),
                         c("blue", "#F9F9F9FF", "red")),
        top_annotation = arm_ha,
        bottom_annotation = column_ha,
        cluster_rows = TRUE, cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        column_names_gp = grid::gpar(fontsize = 8),
        row_names_gp = grid::gpar(fontsize = 8))
dev.off()


#Heatmap of VEH vs TAD for all cell type
################################################
signaling1v2 <- readRDS("signaling1v2.RDS")
signaling1v2 <- signaling1v2[,order(colnames(signaling1v2))]
signaling1v2 <- signaling1v2[, grepl("^(Caf|Gran|Mac|Endo|Tc)", colnames(signaling1v2))]
signaling1v2 <- signaling1v2[rownames(signaling1v2) %in% c("Ccl7","Il22","Ccl12","Il10","Jag2","Jag1","Il6","Csf1","Ccl3","Ccl4","Entpd1"),]

sobj$cell_treat <- paste(sobj$cell_type, sobj$group_id, sep = "_")
cta <- colnames(signaling1v2)
celltype <- gsub("_.*", "", cta)
treat_anno <- gsub(".*_", "", cta)
pal <- hue_pal()(length(unique(celltype)))
palette_key <- setNames(pal, unique(celltype))

pal_arm <- hue_pal()(length(unique(treat_anno)))
pal_arm_key <- setNames(pal_arm, unique(treat_anno))
column_ha <- HeatmapAnnotation(cell_type = celltype,
                               col = list(cell_type = palette_key))
arm_ha <- HeatmapAnnotation(Treatment = treat_anno,
                            col = list(treatment = pal_arm_key))

scaled1v2 <- as.matrix(t(scale(t(signaling1v2))))
scaled1v2 <- na.omit(scaled1v2)
pdf("signaling1v2_heatmap.pdf")
Heatmap(scaled1v2,
        name = "mean scaled\nexpression",
        col = colorRamp2(c(min(scaled1v2), 0, max(scaled1v2)),
                         c("blue", "#F9F9F9FF", "red")),
        top_annotation = arm_ha,
        bottom_annotation = column_ha,
        cluster_rows = TRUE, cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        column_names_gp = grid::gpar(fontsize = 4),
        row_names_gp = grid::gpar(fontsize = 4))
dev.off()


#Heatmap of VEH vs TAD for 3 cell type
######################################
column_list <- c("Granulocyte_TAD","Granulocyte_VEH","Macrophage_TAD","Macrophage_VEH","Tcell_TAD","Tcell_VEH")
threeCelltype <- signaling1v2[,column_list]
sobj$cell_treat <- paste(sobj$cell_type, sobj$group_id, sep = "_")
cta <- colnames(threeCelltype)
celltype <- gsub("_.*", "", cta)
treat_anno <- gsub(".*_", "", cta)
pal <- hue_pal()(length(unique(celltype)))
palette_key <- setNames(pal, unique(celltype))

pal_arm <- hue_pal()(length(unique(treat_anno)))
pal_arm_key <- setNames(pal_arm, unique(treat_anno))
column_ha <- HeatmapAnnotation(cell_type = celltype,
                               col = list(cell_type = palette_key))
arm_ha <- HeatmapAnnotation(Treatment = treat_anno,
                            col = list(treatment = pal_arm_key))

scaled1v2_3cell <- as.matrix(t(scale(t(threeCelltype))))
scaled1v2_3cell <- na.omit(scaled1v2_3cell)
pdf("signaling1v2_heatmap_3celltype.pdf")
Heatmap(scaled1v2_3cell,
        name = "mean scaled\nexpression",
        col = colorRamp2(c(min(scaled1v2_3cell), 0, max(scaled1v2_3cell)),
                         c("blue", "#F9F9F9FF", "red")),
        top_annotation = arm_ha,
        bottom_annotation = column_ha,
        cluster_rows = TRUE, cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        column_names_gp = grid::gpar(fontsize = 4),
        row_names_gp = grid::gpar(fontsize = 4))
dev.off()


#Heatmap of MLIST vs MLIST.TAD for 3 cell type
##############################################
column_list <- c("Granulocyte_MLIST","Granulocyte_MLIST.TAD","Macrophage_MLIST","Macrophage_MLIST.TAD","Tcell_MLIST","Tcell_MLIST.TAD")
threeCelltype <- signaling4v5[,column_list]
sobj$cell_treat <- paste(sobj$cell_type, sobj$group_id, sep = "_")
cta <- colnames(threeCelltype)
celltype <- gsub("_.*", "", cta)
treat_anno <- gsub(".*_", "", cta)
pal <- hue_pal()(length(unique(celltype)))
palette_key <- setNames(pal, unique(celltype))

pal_arm <- hue_pal()(length(unique(treat_anno)))
pal_arm_key <- setNames(pal_arm, unique(treat_anno))
column_ha <- HeatmapAnnotation(cell_type = celltype,
                               col = list(cell_type = palette_key))
arm_ha <- HeatmapAnnotation(Treatment = treat_anno,
                            col = list(treatment = pal_arm_key))

scaled4V5_3cell <- as.matrix(t(scale(t(threeCelltype))))
scaled4V5_3cell <- na.omit(scaled4V5_3cell)
pdf("signaling4V5_heatmap_3celltype.pdf")
Heatmap(scaled4V5_3cell,
        name = "mean scaled\nexpression",
        col = colorRamp2(c(min(scaled4V5_3cell), 0, max(scaled4V5_3cell)),
                         c("blue", "#F9F9F9FF", "red")),
        top_annotation = arm_ha,
        bottom_annotation = column_ha,
        cluster_rows = TRUE, cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        column_names_gp = grid::gpar(fontsize = 4),
        row_names_gp = grid::gpar(fontsize = 4))
dev.off()


##############################################################################
###Heatmap for ligands
Ligang1v2 <- readRDS("Ligand1v2.RDS")
Ligang4v5 <- readRDS("Ligand4v5.RDS")
tadmlisttad <- readRDS("E:/WJHLab/domino/meeting_2/tadmlisttad.RDS")
tadmlisttad <- tadmlisttad[,order(colnames(tadmlisttad))]