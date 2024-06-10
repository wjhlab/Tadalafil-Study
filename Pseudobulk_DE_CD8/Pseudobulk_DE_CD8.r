sessionInfo()

library(scater)

library(Seurat)

library(tidyverse)
library(ggrepel)
library(ggplot2)

library(cowplot)
library(Matrix.utils)

library(edgeR)

library(dplyr)
library(magrittr)

library(Matrix)
library(purrr)
library(reshape2)

library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)

library(RColorBrewer)
library(msigdbr)
library(fgsea)
library('readxl')
library('DESeq2')
#library('ComplexHeatmap')

#library('EnhancedVolcano')
#library('pca3d')
library('sva')

library('DT')

library('fgsea')
library('GSVA')
library(ggplot2)
#library(ggpubr)

library(tidyverse)
library(dplyr)
library(tidyr)
library(data.table)
set.seed(1234)


library(Seurat)
library(RColorBrewer)
library(SeuratData)
library(dplyr)
library(Seurat)
library(patchwork)
library(MAST)
library(org.Hs.eg.db)

#sobj_old <- readRDS("sobj_filtered_transformed_clustered.RDS")
#sample_id = list("G1mus1","G1mus2","G1mus3","G2mus1","G2mus2","G2mus3","G2mus3_rep","G3mus1","G3mus2","G3mus3","G4mus1","G4mus2","G4mus3","G5mus1","G5mus2","G5mus3")
#tube_id = list('1.1','1.2','1.3','2.1','2.2','2.3','2.3.1','3.1','3.2','3.3','4.1','4.2','4.3','5.1','5.2','5.3')
#group = list('VEH','VEH','VEH','TAD','TAD','TAD','TAD','CLIST','CLIST','CLIST','MLIST','MLIST','MLIST','MLIST_TAD','MLIST_TAD','MLIST_TAD')
#index_id = list('B1','B2','B3','B4','B5','B6','C4','B7','B8','B9','B10','B11','B12','C1','C2','C3')

#sobj_old@meta.data['group_id']=NA
#sobj_old@meta.data['sample_id']=NA
#sobj_old@meta.data['tube_id']=NA
#sobj_old@meta.data['index']=NA

#for (i in 1:nrow(sobj_old@meta.data)){
#    indicator = rownames(sobj_old@meta.data)[i]
#    indicator = gsub(".*-","",indicator)
#    sobj_old@meta.data[i,7]=group[as.numeric(indicator)]
#    sobj_old@meta.data[i,8]=sample_id[as.numeric(indicator)]
#    sobj_old@meta.data[i,9]=tube_id[as.numeric(indicator)]
#    sobj_old@meta.data[i,10]=index_id[as.numeric(indicator)]
#}

# Load Data
sobj <- readRDS("filtered_sobj_withSampleid.RDS")

unique(sobj@meta.data$sample_id)

marker <- FindMarkers(sobj, ident.1 = 'Endothelial')

write.csv(marker,"Endothelial_allCells.csv")

marker

sobj$QC <- "mouse_TILs"

sobj@meta.data[1,]
unique(sobj@meta.data$orig.ident)

VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = "QC",pt.size = 0, ncol = 3)

VlnPlot(subset(x=sobj,subset=sample_id=="G5mus1"), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = "QC",pt.size = 0, ncol = 3)

unique(sobj@meta.data$sample_id)

plot1 = FeaturePlot(sobj,features = c("Arg1"), repel = TRUE)
plot2 = FeaturePlot(sobj,features = c("Arg1"), repel = TRUE)
plot3 = FeaturePlot(sobj,features = c("Arg1"), repel = TRUE)
plot4 = FeaturePlot(sobj,features = c("Arg1"), repel = TRUE)
plot5 = FeaturePlot(sobj,features = c("Arg1"), repel = TRUE)
plot6 = FeaturePlot(sobj,features = c("Arg1"), repel = TRUE)

sobj_tcell = subset(x = sobj, subset = cell_type == "Tcell")

sobj_tcell_G45 = subset(x=sobj_tcell,subset = sample_id == c("G4mus1","G4mus2","G4mus3","G5mus1","G5mus2","G5mus3"))



plot4=VlnPlot(sobj_tcell_G45, features = 'Ifng', group.by = "sample_id",slot="scale.data",pt.size = 0)+labs(subtitle='Ifng in G4 G5 Tc')+stat_summary(fun.y = mean, geom='crossbar', colour = "BLACK")
print(plot4)

genes2 = FindMarkers(object=sobj_tcell, group.by = "group_id",ident.1="MLIST",ident.2="MLIST.TAD", test.use = "MAST")

write.csv(genes2,"tcell_2gene.csv")

VEH <- readRDS("veh_2.RDS")
TAD <- readRDS("tad_2.RDS")
MLIST <- readRDS("mlist_2.RDS")
MLISTTAD <- readRDS("mlisttad_2.RDS")

new_tadmlisttad <- readRDS('new_tadmlisttad.rds')
new_tadmlisttad



new_tadmlisttad <- readRDS('new_tadmlisttad.rds')

celltype <- c(levels(sobj))
Ligand <- c()
treat <- list(MLIST=new_tadmlisttad)
cell_treat <- unique(sobj@meta.data$cell_treat)
Ligand4v5 <- data.frame(ligand=NULL,p_val=NULL,avg_log2FC=NULL,pct.1=NULL,pct.2=NULL,p_val_adj=NULL,cell_type=NULL,treatment=NULL,comp=NULL)
for (df_name in names(treat)){
    df <- treat[[df_name]]
    for (i in 1:8){
        temp <- FindMarkers(object = sobj, ident.1 = paste(celltype[i],"TAD",sep="_"), ident.2 = paste(celltype[i],"MLIST.TAD",sep="_"), features = intersect(rownames(sobj), unique(df$ligand)), group.by = "cell_treat",min.pct = 0.0001,min.cells.feature = 0, logfc.threshold = -Inf)
        temp$cell_type <- rep(celltype[i],nrow(temp))
        temp$treatment <- rep(df_name,nrow(temp))
        Ligand <- c(Ligand,rownames(temp))
        Ligand4v5 <- rbind(Ligand4v5,temp)
    }
}
Ligand4v5$ligand <- Ligand

new_list4v5 <- sobj@assays$RNA@data[intersect(rownames(Ligand4v5), rownames(sobj@assays$RNA@data)),]

cta <- unique(sobj$cell_treat)[grepl("TAD$|MLIST.TAD$", unique(sobj$cell_treat))]
signaling4v5 <- matrix(0, ncol = length(cta),
                      nrow = length(rownames(new_list4v5)))
rownames(signaling4v5) <- rownames(new_list4v5)
colnames(signaling4v5) <- cta
for(clust in cta){
    sig <- rowMeans(new_list4v5[, which(sobj$cell_treat == clust)])
    # ensure the rowmeans order matches the rownames of the signaling matrix
    sig <- sig[rownames(signaling4v5)]
    signaling4v5[, clust] <- sig
}

saveRDS(signaling4v5,file = 'tadmlisttad.RDS')

mlist1mlisttad0 <- readRDS('new_mlist1mlisttad0.rds')
mlisttad1mlist0 <- readRDS('new_mlisttad1mlist0.rds')

celltype <- c(levels(sobj))
Ligand <- c()
treat <- list(MLIST=mlist1mlisttad0,MLISTTAD=mlisttad1mlist0)
cell_treat <- unique(sobj@meta.data$cell_treat)
Ligand4v5 <- data.frame(ligand=NULL,p_val=NULL,avg_log2FC=NULL,pct.1=NULL,pct.2=NULL,p_val_adj=NULL,cell_type=NULL,treatment=NULL,comp=NULL)
for (df_name in names(treat)){
    df <- treat[[df_name]]
    for (i in 1:8){
        temp <- FindMarkers(object = sobj, ident.1 = paste(celltype[i],"MLIST",sep="_"), ident.2 = paste(celltype[i],"MLIST.TAD",sep="_"), features = intersect(rownames(sobj), unique(df$ligand)), group.by = "cell_treat",min.pct = 0.0001,min.cells.feature = 0, logfc.threshold = -Inf)
        temp$cell_type <- rep(celltype[i],nrow(temp))
        temp$treatment <- rep(df_name,nrow(temp))
        Ligand <- c(Ligand,rownames(temp))
        Ligand4v5 <- rbind(Ligand4v5,temp)
    }
}
Ligand4v5$ligand <- Ligand

new_list4v5 <- sobj@assays$RNA@data[intersect(rownames(Ligand4v5), rownames(sobj@assays$RNA@data)),]

cta <- unique(sobj$cell_treat)[grepl("MLIST$|MLIST.TAD$", unique(sobj$cell_treat))]
signaling4v5 <- matrix(0, ncol = length(cta),
                      nrow = length(rownames(new_list4v5)))
rownames(signaling4v5) <- rownames(new_list4v5)
colnames(signaling4v5) <- cta
for(clust in cta){
    sig <- rowMeans(new_list4v5[, which(sobj$cell_treat == clust)])
    # ensure the rowmeans order matches the rownames of the signaling matrix
    sig <- sig[rownames(signaling4v5)]
    signaling4v5[, clust] <- sig
}

saveRDS(signaling4v5,file = 'signaling4v5.RDS')





sobj$cell_treat <- paste(sobj$cell_type, sobj$group_id, sep = "_")

unique(sobj$cell_treat)

### All included ligand list

celltype <- c(levels(sobj))
Ligand <- c()
treat <- list(VEH=VEH,TAD=TAD)
cell_treat <- unique(sobj@meta.data$cell_treat)
Ligand1v2 <- data.frame(ligand=NULL,p_val=NULL,avg_log2FC=NULL,pct.1=NULL,pct.2=NULL,p_val_adj=NULL,cell_type=NULL,treatment=NULL,comp=NULL)
for (df_name in names(treat)){
    df <- treat[[df_name]]
    for (i in 1:8){
        temp <- FindMarkers(object = sobj, ident.1 = paste(celltype[i],"VEH",sep="_"), ident.2 = paste(celltype[i],"TAD",sep="_"), features = intersect(rownames(sobj), unique(df$ligand)), group.by = "cell_treat",min.pct = 0.0001,min.cells.feature = 0, logfc.threshold = -Inf)
        temp$cell_type <- rep(celltype[i],nrow(temp))
        temp$treatment <- rep(df_name,nrow(temp))
        Ligand <- c(Ligand,rownames(temp))
        Ligand1v2 <- rbind(Ligand1v2,temp)
    }
}
Ligand1v2$ligand <- Ligand

celltype <- c(levels(sobj))
Ligand <- c()
treat <- list(MLIST=MLIST,MLISTTAD=MLISTTAD)
cell_treat <- unique(sobj@meta.data$cell_treat)
Ligand4v5 <- data.frame(ligand=NULL,p_val=NULL,avg_log2FC=NULL,pct.1=NULL,pct.2=NULL,p_val_adj=NULL,cell_type=NULL,treatment=NULL,comp=NULL)
for (df_name in names(treat)){
    df <- treat[[df_name]]
    for (i in 1:8){
        temp <- FindMarkers(object = sobj, ident.1 = paste(celltype[i],"MLIST",sep="_"), ident.2 = paste(celltype[i],"MLIST.TAD",sep="_"), features = intersect(rownames(sobj), unique(df$ligand)), group.by = "cell_treat",min.pct = 0.0001,min.cells.feature = 0, logfc.threshold = -Inf)
        temp$cell_type <- rep(celltype[i],nrow(temp))
        temp$treatment <- rep(df_name,nrow(temp))
        Ligand <- c(Ligand,rownames(temp))

        Ligand4v5 <- rbind(Ligand4v5,temp)
    }
}
Ligand4v5$ligand <- Ligand

dim(Ligand1v2)
dim(Ligand4v5)

sigLigand1v2 <- Ligand1v2[abs(Ligand1v2$avg_log2FC)>0.5 & Ligand1v2$p_val_adj < 0.05,]
sigLigand4v5 <- Ligand4v5[abs(Ligand4v5$avg_log2FC)>0.5 & Ligand4v5$p_val_adj < 0.05,]
Ligandlist <- c(rownames(Ligand1v2),rownames(Ligand4v5))
dim(Ligand1v2)
dim(Ligand4v5)
length(Ligandlist)

new_list1v2 <- sobj@assays$RNA@data[intersect(rownames(Ligand1v2), rownames(sobj@assays$RNA@data)),]
new_list4v5 <- sobj@assays$RNA@data[intersect(rownames(Ligand4v5), rownames(sobj@assays$RNA@data)),]

#cta <- unique(sobj$cell_treat)[order(unique(sobj$cell_treat))]
cta <- unique(sobj$cell_treat)[grepl("_[VT]AD$|_[VT]EH$", unique(sobj$cell_treat))]
signaling1v2 <- matrix(0, ncol = length(cta),
                      nrow = length(rownames(new_list1v2)))
rownames(signaling1v2) <- rownames(new_list1v2)
colnames(signaling1v2) <- cta
for(clust in cta){
    sig <- rowMeans(new_list1v2[, which(sobj$cell_treat == clust)])
    # ensure the rowmeans order matches the rownames of the signaling matrix
    sig <- sig[rownames(signaling1v2)]
    signaling1v2[, clust] <- sig
}

cta <- unique(sobj$cell_treat)[grepl("MLIST", unique(sobj$cell_treat))]
signaling4v5 <- matrix(0, ncol = length(cta),
                      nrow = length(rownames(new_list4v5)))
rownames(signaling4v5) <- rownames(new_list4v5)
colnames(signaling4v5) <- cta
for(clust in cta){
    sig <- rowMeans(new_list4v5[, which(sobj$cell_treat == clust)])
    # ensure the rowmeans order matches the rownames of the signaling matrix
    sig <- sig[rownames(signaling4v5)]
    signaling4v5[, clust] <- sig
}

"csf1" %in% rownames(signaling1v2)

cta <- unique(sobj$cell_treat)[grepl("_[VT]AD$|_[VT]EH$", unique(sobj$cell_treat))]
cta

saveRDS(signaling1v2,"signaling1v2.RDS")
saveRDS(signaling4v5,"signaling4v5.RDS")

saveRDS(signaling1v2,"signaling1v2.RDS")
saveRDS(signaling4v5,"signaling4v5.RDS")
saveRDS(sigLigand1v2,"sigLigand1v2.RDS")
saveRDS(sigLigand4v5,"sigLigand4v5.RDS")
saveRDS(Ligand1v2,"Ligand1v2.RDS")
saveRDS(Ligand4v5,"Ligand4v5.RDS")
saveRDS(new_list1v2,"new_list1v2.RDS")
saveRDS(new_list4v5,"new_list4v5.RDS")

intersect(Ligandlist,colnames(sobj@assays$RNA@data))

### Applying Jacob's code to get the percent of expression of ligand by receptor
cl_rec_percent_mlist = NULL
ser_receptors = unique(mlistTAD_2$receptor)
for(rec in ser_receptors){
  rec_percent <- sapply(
    X = unique(Idents(sobj)),
    FUN = function(x){
      # percentage of cells in cluster with non-zero expression of receptor
      sum(dom_mlist@counts[rec,dom_mlist@clusters == x] > 0) / length(dom_mlist@counts[rec,dom_mlist@clusters == x])
    }
  )
  cl_rec_percent_mlist <- rbind(cl_rec_percent_mlist, rec_percent)
}
rownames(cl_rec_percent_mlist) = ser_receptors

cell_clusters <- Idents(sobj)
cluster_names <- names(table(cell_clusters))
receptor <- vehTAD_2$recptor

# Get the expression matrix from the 'RNA' assay
gene_expression <- sobj@assays$RNA@data

# Get the index of the receptor gene in the expression matrix
receptor_index <- which(rownames(gene_expression) == receptor)

# Create an empty data frame to store the expression values
gene_expression_table <- data.frame(Cluster = character(),
                                    Gene = character(),
                                    Expression = logical(),
                                    stringsAsFactors = FALSE)

for (cluster_id in unique(cell_clusters)) {
  cluster_expression <- gene_expression[receptor_index, cell_clusters == cluster_id]
  is_expressed <- any(cluster_expression > 0)
  cluster_name <- cluster_names[cluster_id]
  
  # Add the expression values to the table
  gene_expression_table <- rbind(gene_expression_table,
                                 data.frame(Cluster = cluster_name,
                                            Gene = receptor,
                                            Expression = is_expressed,
                                            stringsAsFactors = FALSE))
}

#gene_expression <- sobj@assays$RNA@data
#trem2_expression <- gene_expression["Trem2", ]
#VlnPlot(object = sobj, features = 'Trem2')

cell_clusters <- Idents(sobj)

for (cluster_id in unique(cell_clusters)) {
  cluster_expression <- gene_expression["Trem2", cell_clusters == cluster_id]
  is_expressed <- any(cluster_expression > 0)
  print(paste("Cluster", cluster_id, "Trem2 expression:", is_expressed))
}
sum(dom@counts[rec,dom@clusters == x] > 0) / length(dom@counts[rec,dom@clusters == x])

cl_rec_percent = NULL
ser_receptors = vehTAD_2$recptor
for(rec in ser_receptors){
        rec_percent <- sapply(
          X = Idents(sobj),
          FUN = function(x){
            # percentage of cells in cluster with non-zero expression of receptor
            sum(dom@counts[rec,dom@clusters == x] > 0) / length(dom@counts[rec,dom@clusters == x])
          }
        )
        cl_rec_percent <- rbind(cl_rec_percent, rec_percent)
      }

#saveRDS(sobj, "filtered_sobj_withSampleid.RDS")

colnames(sobj@meta.data)

unique(sobj@meta.data$cell_type)

no_tumor = subset(x = sobj, subset = cell_type != c("Tumor"))

pdf("updated_umap.pdf")
DimPlot(no_tumor, group.by = "cell_type",label = TRUE) + NoLegend()
dev.off()

##Needed?
#avg_exp <- AverageExpression(seur, group.by = "Cell_Type", slot = "counts") # this may be used for downstream logFC
#exp <- avg_exp$Spatial %>% as_tibble(rownames = "symbol") %>% dplyr::select(symbol, DC) %>% rename(mean_exp = DC)

## Load Pathways
##Get Pathways
hallmark_df = msigdbr(species = "Mus musculus", category = "H")
hallmark_list = hallmark_df %>% split(x = .$gene_symbol, f = .$gs_name)

oncogenic_df = msigdbr(species = "Mus musculus", category = "C6")
oncogenic_list = oncogenic_df %>% split(x = .$gene_symbol, f = .$gs_name)

immunologic_df = msigdbr(species = "Mus musculus", category = "C7")
immunologic_list = immunologic_df %>% split(x = .$gene_symbol, f = .$gs_name)

Curated_df = msigdbr(species = "Mus musculus", category = "C2")
kegg_list = Curated_df[Curated_df$gs_subcat=="CP:KEGG",] %>% split(x = .$gene_symbol, f = .$gs_name)
reactome_list = Curated_df[Curated_df$gs_subcat=="CP:REACTOME",] %>% split(x = .$gene_symbol, f = .$gs_name)
NABA_list = Curated_df[Curated_df$gs_name%in%c("NABA_BASEMENT_MEMBRANES","NABA_COLLAGENS","NABA_CORE_MATRISOME","NABA_ECM_AFFILIATED","NABA_ECM_GLYCOPROTEINS","NABA_ECM_REGULATORS","NABA_MATRISOME","NABA_MATRISOME_ASSOCIATED","NABA_PROTEOGLYCANS","NABA_SECRETED_FACTORS"),] %>% split(x = .$gene_symbol, f = .$gs_name)

Ontology_df = msigdbr(species = "Mus musculus", category = "C5")
GOBP_list = Ontology_df[Ontology_df$gs_subcat=="GO:BP",] %>% split(x = .$gene_symbol, f = .$gs_name)
GOMF_list = Ontology_df[Ontology_df$gs_subcat=="GO:MF",] %>% split(x = .$gene_symbol, f = .$gs_name)


groups <- sobj@meta.data[, c("sample_id",  "cell_type", "group_id")]

unique(groups$sample_id)

pb <- aggregate.Matrix(t(sobj@assays$RNA@counts), # raw RNA counts
                       groupings = groups, fun = "sum")

pb_mtx <- t(as.matrix(pb))

splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 6), 
                 `[`, 2)

pb_list <- split.data.frame(pb, 
                       factor(splitf)) %>%
        lapply(function(u) 
                set_colnames(t(u), 
                             stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

groups <- sobj@meta.data[, c("sample_id",  "cell_type", "group_id")]
pb <- aggregate.Matrix(t(sobj@assays$RNA@counts), # raw RNA counts
                       groupings = groups, fun = "sum")



sampAnnot <- data.frame(sample_id = rownames(pb))
rownames(sampAnnot) <- sampAnnot$sample_id


sampAnnot$sample_id <- data.frame(str_split_fixed(rownames(sampAnnot), "_", 6))$X1    #id
sampAnnot$cell_type <- data.frame(str_split_fixed(rownames(sampAnnot), "_", 6))$X2    #cellType
sampAnnot$group_id <- data.frame(str_split_fixed(rownames(sampAnnot), "_", 6))$X3


pb_mtx <- t(as.matrix(pb))
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 6), 
                 `[`, 2)
pb_list <- split.data.frame(pb, 
                       factor(splitf)) %>%
        lapply(function(u) 
                set_colnames(t(u), 
                             stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

#saveRDS(pb, "sobj_pseudobulk.rds")

pb

## Evaluation of sample quality based upon distribution of reads
boxplot(log2(pb_mtx+1), ylab='log2(counts+1)', 
        xlab='samples', las=2, cex.axis = 0.25, 
        border = ifelse(apply(log2(pb_mtx+1),2,quantile,probs=0.75)>0,
                        'black','red'))
legend('topleft', pch=c('-'), legend = c('upperquant > 0', 'upperquant < 0'),
       col=c('black','red'), cex=0.5)

datatable(cbind(sampAnnot,
                lowCountFilter=apply(log2(pb_mtx+1),2,quantile,probs=0.75)>0))

## Subset based on upperquant >0
sampAnnotQuantFilter <- sampAnnot[apply(log2(pb_mtx+1),2,quantile,probs=0.75)>0,]
pb_mtx_quant <- pb_mtx[,apply(log2(pb_mtx+1),2,quantile,probs=0.75)>0]



#[`DESeq2`](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) is used for variance stabilization transformation for count normalization. 

## Obtain normalized log counts with vst 
dds <- DESeqDataSetFromMatrix(countData = round(pb_mtx_quant),
                              colData = sampAnnotQuantFilter,
                              design=~group_id)
vstDat <- assay(vst(dds))
vstDat <- vstDat[apply(vstDat,1,sd)>0,]

## Perform PCA analysis to evaluate sample clustering

#We use PCA of the `vst` normalized data to evaluate sample clustering by cell type. We observe largely robust clustering by Batch.

pcs <- prcomp(t(vstDat),scale=T)
pca2d(pcs,
      group=sampAnnotQuantFilter$Batch, legend="topleft")
title('PCA of vst counts after quantile filtering')

pca2d(pcs,
      group=sampAnnotQuantFilter$CellType, legend="topleft")
title('PCA of vst counts after quantile filtering')

# Evaluation of batch effects
pca2d(pcs,
      group=sampAnnotQuantFilter$Batch, legend="topleft")
title('PCA of vst counts of by batch')

datCombat <- sva::ComBat(vstDat,
                         batch=sampAnnotQuantFilter$Batch,
                         mod=model.matrix(~sampAnnotQuantFilter$Arm))

pcsCombat <- prcomp(t(datCombat),scale=T)
pca2d(pcsCombat,
      group=sampAnnotQuantFilter$CellType, legend="topright")
title('PCA of vst counts of batch corrected data')
pca2d(pcsCombat,
      group=sampAnnotQuantFilter$Batch, legend="topright")
title('PCA of vst counts of batch corrected data')

# Final CellType Numbers
table(sampAnnotQuantFilter$CellType, sampAnnotQuantFilter$Arm)
table(sampAnnotQuantFilter$CellType, sampAnnotQuantFilter$OS)
table(sampAnnotQuantFilter$OS, sampAnnotQuantFilter$Arm)
table(sampAnnotQuantFilter$Arm)
table(sampAnnotQuantFilter$OS)





### DE analysis

# Set up pseudobulk and DE function
ARM_DE <- function(CELLTYPE){

##Subset to celltype of interest
meta <- sampAnnotQuantFilter[sampAnnotQuantFilter$cell_type==CELLTYPE,]
counts <- pb_mtx_quant[,rownames(meta)]

##Make dds
dds <- DESeqDataSetFromMatrix(round(counts), 
                              colData = meta, 
                              design = ~group_id)

dds <- DESeq(dds)

##Make Results Comparisons
one_twoRES <- results(dds, contrast = c("group_id","VEH","TAD"))
four_fiveRES <- results(dds, contrast = c("group_id",'MLIST','MLIST.TAD'))


##LfcShrink -- ashr
#one_twoRES <- lfcShrink(dds,contrast = c("group_id","VEH","TAD"), 
#                   res=one_twoRES, type="ashr")
#four_fiveRES <- lfcShrink(dds,contrast = c("group_id",'MLIST','MLIST.TAD'), 
#                   res=four_fiveRES, type="ashr")


##Filter out NAs
one_twoRES <- one_twoRES[!is.na(one_twoRES$padj),]
four_fiveRES <- four_fiveRES[!is.na(four_fiveRES$padj),]


##Extract Stats
one_twoSTATS <- one_twoRES$log2FoldChange
names(one_twoSTATS) <- rownames(one_twoRES)

four_fiveSTATS <- four_fiveRES$log2FoldChange
names(four_fiveSTATS) <- rownames(four_fiveRES)

##Run Pathway Analysis
ONEvTWOgsResults <- list(HALLMARK=fgsea(pathways=hallmark_list, 
                               stats=one_twoSTATS),
                    KEGG=fgsea(pathways=kegg_list, 
                               stats=one_twoSTATS),
                    REACTOME=fgsea(pathways=reactome_list,
                                   stats=one_twoSTATS), 
                    IMMUNESIG=fgsea(pathways=immunologic_list,
                                     stats=one_twoSTATS),
                    GOBP = fgsea(pathways=GOBP_list,
                                     stats=one_twoSTATS),
                    GOMF = fgsea(pathways=GOMF_list,
                                 stats=one_twoSTATS),
                    NABA = fgsea(pathways=NABA_list,
                                 stats=one_twoSTATS))

FOURvFIVEgsResults <- list(HALLMARK=fgsea(pathways=hallmark_list, 
                               stats=four_fiveSTATS),
                    KEGG=fgsea(pathways=kegg_list, 
                               stats=four_fiveSTATS),
                    REACTOME=fgsea(pathways=reactome_list,
                                   stats=four_fiveSTATS), 
                    IMMUNESIG=fgsea(pathways=immunologic_list,
                                     stats=four_fiveSTATS),
                    GOBP = fgsea(pathways=GOBP_list,
                                     stats=four_fiveSTATS),
                    GOMF = fgsea(pathways=GOMF_list,
                                 stats=four_fiveSTATS),
                    NABA = fgsea(pathways=NABA_list,
                                 stats=four_fiveSTATS))

vst_counts = vst(dds)  

  return(list(one_twoRES=one_twoRES, ONEvTWOgs=ONEvTWOgsResults, four_fiveRES=four_fiveRES, FOURvFIVEgs=FOURvFIVEgsResults, Meta = meta, vst_counts=vst_counts))
}




## Macrophage
Mac_DE <- ARM_DE("Macrophage")

##1v2
Mac_1v2DE <- Mac_DE$one_twoRES

Mac_ONEvTWO_DEGenes <- row.names(Mac_1v2DE)[Mac_1v2DE$padj < 0.05 &
                              abs(Mac_1v2DE$log2FoldChange) > 0.5]

Mac_ONEvTWO_TopDE <- Mac_1v2DE[Mac_1v2DE$padj < 0.05 &
                              abs(Mac_1v2DE$log2FoldChange) > 0.5,]

write.csv(Mac_ONEvTWO_DEGenes, "Macrophage_TopDE_1v2filtered.csv")
write.csv(Mac_1v2DE,"Macrophage_TopDE_1v2unfiltered.csv")

##4v5
Mac_4v5DE <- Mac_DE$four_fiveRES

Mac_FOURvFIVE_DEGenes <- row.names(Mac_4v5DE)[Mac_4v5DE$padj < 0.05 &
                              abs(Mac_4v5DE$log2FoldChange) > 0.5]

Mac_FOURvFIVE_TopDE <- Mac_4v5DE[Mac_4v5DE$padj < 0.05 &
                              abs(Mac_4v5DE$log2FoldChange) > 0.5,]

write.csv(Mac_FOURvFIVE_TopDE, "Macrophage_TopDE_4v5filtered.csv")
write.csv(Mac_4v5DE, "Macrophage_TopDE_4v5unfiltered.csv")

##Extract Metadata
meta_Mac <- Mac_DE$Meta
ONEvTWOmeta_Mac <- meta_Mac[!meta_Mac$group_id=="CLIST"&
                            !meta_Mac$group_id=="MLIST"&
                            !meta_Mac$group_id=="MLIST.TAD",]
FOURvFIVEmeta_Mac <- meta_Mac[!meta_Mac$group_id=="CLIST"&
                              !meta_Mac$group_id=="VEH"&
                              !meta_Mac$group_id=="TAD",]

## Granulocyte
Gra_DE <- ARM_DE("Granulocyte")

##1v2
Gra_1v2DE <- Gra_DE$one_twoRES

Gra_ONEvTWO_DEGenes <- row.names(Gra_1v2DE)[Gra_1v2DE$padj < 0.05 &
                              abs(Gra_1v2DE$log2FoldChange) > 0.5]

Gra_ONEvTWO_TopDE <- Gra_1v2DE[Gra_1v2DE$padj < 0.05 &
                              abs(Gra_1v2DE$log2FoldChange) > 0.5,]

write.csv(Gra_ONEvTWO_DEGenes, "Granulocyte_TopDE_1v2filtered.csv")
write.csv(Gra_1v2DE,"Granulocyte_TopDE_1v2unfiltered.csv")

##4v5
Gra_4v5DE <- Gra_DE$four_fiveRES

Gra_FOURvFIVE_DEGenes <- row.names(Gra_4v5DE)[Gra_4v5DE$padj < 0.05 &
                              abs(Gra_4v5DE$log2FoldChange) > 0.5]

Gra_FOURvFIVE_TopDE <- Gra_4v5DE[Gra_4v5DE$padj < 0.05 &
                              abs(Gra_4v5DE$log2FoldChange) > 0.5,]

write.csv(Gra_FOURvFIVE_TopDE, "Granulocyte_TopDE_4v5filtered.csv")
write.csv(Gra_4v5DE, "Granulocyte_TopDE_4v5unfiltered.csv")

##Extract Metadata
meta_Gra <- Gra_DE$Meta
ONEvTWOmeta_Gra <- meta_Gra[!meta_Gra$group_id=="CLIST"&
                            !meta_Gra$group_id=="MLIST"&
                            !meta_Gra$group_id=="MLIST.TAD",]
FOURvFIVEmeta_Gra <- meta_Gra[!meta_Gra$group_id=="CLIST"&
                              !meta_Gra$group_id=="VEH"&
                              !meta_Gra$group_id=="TAD",]


## Tcell
Tcell_DE <- ARM_DE("Tcell")

##1v2
Tcell_1v2DE <- Tcell_DE$one_twoRES

Tcell_ONEvTWO_DEGenes <- row.names(Tcell_1v2DE)[Tcell_1v2DE$padj < 0.05 &
                              abs(Tcell_1v2DE$log2FoldChange) > 0.5]

Tcell_ONEvTWO_TopDE <- Tcell_1v2DE[Tcell_1v2DE$padj < 0.05 &
                              abs(Tcell_1v2DE$log2FoldChange) > 0.5,]

write.csv(Tcell_ONEvTWO_DEGenes, "Tcell_TopDE_1v2filtered.csv")
write.csv(Tcell_1v2DE,"Tcell_TopDE_1v2unfiltered.csv")

##4v5
Tcell_4v5DE <- Tcell_DE$four_fiveRES

Tcell_FOURvFIVE_DEGenes <- row.names(Tcell_4v5DE)[Tcell_4v5DE$padj < 0.05 &
                              abs(Tcell_4v5DE$log2FoldChange) > 0.5]

Tcell_FOURvFIVE_TopDE <- Tcell_4v5DE[Tcell_4v5DE$padj < 0.05 &
                              abs(Tcell_4v5DE$log2FoldChange) > 0.5,]

write.csv(Tcell_FOURvFIVE_TopDE, "Tcell_TopDE_4v5filtered.csv")
write.csv(Tcell_4v5DE, "Tcell_TopDE_4v5unfiltered.csv")

##Extract Metadata
meta_Tcell <- Tcell_DE$Meta
ONEvTWOmeta_Tcell <- meta_Tcell[!meta_Tcell$group_id=="CLIST"&
                            !meta_Tcell$group_id=="MLIST"&
                            !meta_Tcell$group_id=="MLIST.TAD",]
FOURvFIVEmeta_Tcell <- meta_Tcell[!meta_Tcell$group_id=="CLIST"&
                              !meta_Tcell$group_id=="VEH"&
                              !meta_Tcell$group_id=="TAD",]



df = data.frame(pathway=NULL,pval=NULL,padj=NULL,log2err=NULL,ES=NULL,NES=NULL,size=NULL,leadingEdge=NULL)
for (i in 1:length(Tcell_DE$ONEvTWOgs)){
    a = data.frame(Tcell_DE$ONEvTWOgs[i])
    colnames(a) = c("pathway","pval","padj","log2err","ES","NES","size","leadingEdge")
    df = rbind(df,a)
}
df$leadingEdge = as.vector(df$leadingEdge,'character')
write.csv(df,file = "Tcell_1v2.csv")

df = data.frame(pathway=NULL,pval=NULL,padj=NULL,log2err=NULL,ES=NULL,NES=NULL,size=NULL,leadingEdge=NULL)
for (i in 1:length(Tcell_DE$FOURvFIVEgs)){
    a = data.frame(Tcell_DE$FOURvFIVEgs[i])
    colnames(a) = c("pathway","pval","padj","log2err","ES","NES","size","leadingEdge")
    df = rbind(df,a)
}
df$leadingEdge = as.vector(df$leadingEdge,'character')
write.csv(df,file = "Tcell_4v5.csv")

df = data.frame(pathway=NULL,pval=NULL,padj=NULL,log2err=NULL,ES=NULL,NES=NULL,size=NULL,leadingEdge=NULL)
for (i in 1:length(Gra_DE$ONEvTWOgs)){
    a = data.frame(Gra_DE$ONEvTWOgs[i])
    colnames(a) = c("pathway","pval","padj","log2err","ES","NES","size","leadingEdge")
    df = rbind(df,a)
}
df$leadingEdge = as.vector(df$leadingEdge,'character')
write.csv(df,file = "Gra_1v2.csv")

df = data.frame(pathway=NULL,pval=NULL,padj=NULL,log2err=NULL,ES=NULL,NES=NULL,size=NULL,leadingEdge=NULL)
for (i in 1:length(Gra_DE$FOURvFIVEgs)){
    a = data.frame(Gra_DE$FOURvFIVEgs[i])
    colnames(a) = c("pathway","pval","padj","log2err","ES","NES","size","leadingEdge")
    df = rbind(df,a)
}
df$leadingEdge = as.vector(df$leadingEdge,'character')
write.csv(df,file = "Gra_4v5.csv")

df = data.frame(pathway=NULL,pval=NULL,padj=NULL,log2err=NULL,ES=NULL,NES=NULL,size=NULL,leadingEdge=NULL)
for (i in 1:length(Mac_DE$ONEvTWOgs)){
    a = data.frame(Mac_DE$ONEvTWOgs[i])
    colnames(a) = c("pathway","pval","padj","log2err","ES","NES","size","leadingEdge")
    df = rbind(df,a)
}
df$leadingEdge = as.vector(df$leadingEdge,'character')
write.csv(df,file = "Tcell_1v2.csv")

df = data.frame(pathway=NULL,pval=NULL,padj=NULL,log2err=NULL,ES=NULL,NES=NULL,size=NULL,leadingEdge=NULL)
for (i in 1:length(Mac_DE$FOURvFIVEgs)){
    a = data.frame(Mac_DE$FOURvFIVEgs[i])
    colnames(a) = c("pathway","pval","padj","log2err","ES","NES","size","leadingEdge")
    df = rbind(df,a)
}
df$leadingEdge = as.vector(df$leadingEdge,'character')
write.csv(df,file = "Mac_4v5.csv")





colnames(a) = c("pathway","pval","padj","log2err","ES","NES","size","leadingEdge")

fwrite(as.data.frame(Mac_DE$ONEvTWOgs), file="Mac_1v2.xlsx")
fwrite(as.data.frame(Mac_DE$FOURvFIVEgs), file="Mac_4v5.xlsx")
fwrite(as.data.frame(Gra_DE$ONEvTWOgs), file="Gra_1v2.xlsx")
fwrite(as.data.frame(Gra_DE$FOURvFIVEgs), file="Gra_4v5.xlsx")
fwrite(as.data.frame(Tcell_DE$ONEvTWOgs), file="Tcell_1v2.xlsx")
fwrite(as.data.frame(Tcell_DE$FOURvFIVEgs), file="Tcell_4v5.xlsx")

genelist_gra1v2_hallmark <- list(HALLMARK_INFLAMMATORY_RESPONSE =
                      unlist(Gra_DE$ONEvTWOgs$HALLMARK$leadingEdge[which(Gra_DE$ONEvTWOgs$HALLMARK$pathway == 'HALLMARK_INFLAMMATORY_RESPONSE')]),
                             HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION =
                      unlist(Gra_DE$ONEvTWOgs$HALLMARK$leadingEdge[which(Gra_DE$ONEvTWOgs$HALLMARK$pathway == 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION')]),
                             HALLMARK_IL2_STAT5_SIGNALING =
                      unlist(Gra_DE$ONEvTWOgs$HALLMARK$leadingEdge[which(Gra_DE$ONEvTWOgs$HALLMARK$pathway == 'HALLMARK_IL2_STAT5_SIGNALING')]))

sobj_gra1v2_hallmark <- AddModuleScore(sobj, genelist_gra1v2_hallmark, pool = NULL, nbin = 24, ctrl = 100,
                       k = FALSE, assay = NULL, name = "module", seed = 1)

colnames(sobj_gra1v2_hallmark@meta.data)[grepl("module",colnames(sobj_gra1v2_hallmark@meta.data))] <- names(genelist_gra1v2_hallmark)

head(sobj_gra1v2_hallmark@meta.data)

module_matrix = as.matrix(sobj_gra1v2_hallmark@meta.data)
module_matrix = module_matrix[,c('HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','HALLMARK_IL2_STAT5_SIGNALING')]
matrix_num <- matrix(as.numeric(module_matrix),    # Convert to numeric matrix
                  ncol = ncol(module_matrix))
row.names(matrix_num) = row.names(module_matrix)
colnames(matrix_num) = colnames(module_matrix)
head(matrix_num)

groups <- sobj_gra1v2_hallmark@meta.data[, c("sample_id",  "cell_type", "group_id")]
pb <- aggregate.Matrix(matrix_num,groupings = groups, fun = "mean")


class(pb)

pheatmap(pb)





genelist_gra4v5_hallmark <- list(HALLMARK_IL2_STAT5_SIGNALING=
                      Gra_DE$FOURvFIVEgs$HALLMARK$leadingEdge[which(Gra_DE$FOURvFIVEgs$HALLMARK$pathway=='HALLMARK_IL2_STAT5_SIGNALING')],
                                 HALLMARK_INFLAMMATORY_RESPONSE=
                      Gra_DE$FOURvFIVEgs$HALLMARK$leadingEdge[which(Gra_DE$FOURvFIVEgs$HALLMARK$pathway=='HALLMARK_INFLAMMATORY_RESPONSE')],
                                 HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY=
                      Gra_DE$FOURvFIVEgs$HALLMARK$leadingEdge[which(Gra_DE$FOURvFIVEgs$HALLMARK$pathway=='HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY')],
                                HALLMARK_APOPTOSIS=
                      Gra_DE$FOURvFIVEgs$HALLMARK$leadingEdge[which(Gra_DE$FOURvFIVEgs$HALLMARK$pathway=='HALLMARK_APOPTOSIS')],
                                HALLMARK_TNFA_SIGNALING_VIA_NFKB=
                      Gra_DE$FOURvFIVEgs$HALLMARK$leadingEdge[which(Gra_DE$FOURvFIVEgs$HALLMARK$pathway=='HALLMARK_TNFA_SIGNALING_VIA_NFKB')],
                                HALLMARK_INTERFERON_ALPHA_RESPONSE=
                      Gra_DE$FOURvFIVEgs$HALLMARK$leadingEdge[which(Gra_DE$FOURvFIVEgs$HALLMARK$pathway=='HALLMARK_INTERFERON_ALPHA_RESPONSE')],
                                HALLMARK_INTERFERON_GAMMA_RESPONSE=
                      Gra_DE$FOURvFIVEgs$HALLMARK$leadingEdge[which(Gra_DE$FOURvFIVEgs$HALLMARK$pathway=='HALLMARK_INTERFERON_GAMMA_RESPONSE')],
                                HALLMARK_HYPOXIA=
                      Gra_DE$FOURvFIVEgs$HALLMARK$leadingEdge[which(Gra_DE$FOURvFIVEgs$HALLMARK$pathway=='HALLMARK_HYPOXIA')])

Gra_DE$FOURvFIVEgs$HALLMARK$leadingEdge[which(Gra_DE$FOURvFIVEgs$HALLMARK$pathway=='HALLMARK_IL2_STAT5_SIGNALING')]

sobj_gra4v5_hallmark <- AddModuleScore(sobj, genelist_gra4v5_hallmark, pool = NULL, nbin = 24, ctrl = 100,
                       k = FALSE, assay = NULL, name = "module", seed = 1)

sobj_gra4v5_hallmark@meta.data

write.csv(Tcell_DE$four_fiveRES,"test.csv")

saveRDS(Gra_DE,"Gra_DE2.rds")
saveRDS(Mac_DE,"Mac_DE2.rds")
saveRDS(Tcell_DE,"Tcell_DE2.rds")



a=data.frame(Mac_4v5DE)
if(a$padj<0.05){
    a$new_category="nontrivial"
} else {
    a$new_category="ns"
}
get_point_color <- function(x) {
  dplyr::case_when(x <= 0.05 ~ "black",
                   x > 0.05 ~ "red")
}
options(repr.plot.width =15, repr.plot.height =15)

pdf("volPlot_mac_4v5.pdf",vol_plot)
vol_plot <- 
  ggplot(data=a,aes(x = log2FoldChange,
             y = -log10(padj),colour = get_point_color(a$padj))) + 
  geom_point(size=2) +
  geom_text(size=6,aes(label=ifelse(padj<0.05,as.character(rownames(a)),'')),hjust=0,vjust=0,) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  NoLegend() +
  ggtitle("vol_plot for Macrophage MLIST vs MLIST.TAD") +
  theme(plot.title = element_text(size = 20, face = "bold"))

vol_plot
dev.off()


a=data.frame(Mac_4v5DE)
if(a$padj<0.05){
    a$new_category="nontrivial"
} else {
    a$new_category="ns"
}
get_point_color <- function(x) {
  dplyr::case_when(x <= 0.05 ~ "black",
                   x > 0.05 ~ "red")
}
options(repr.plot.width =15, repr.plot.height =15)

vol_plot <- 
  ggplot(data=a,aes(x = log2FoldChange,
             y = -log10(padj),colour = get_point_color(a$padj))) + 
  geom_point(size=2) +
  geom_text(size=6,aes(label=ifelse(padj<0.05,as.character(rownames(a)),'')),hjust=0,vjust=0,) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  NoLegend() +
  ggtitle("vol_plot for Macrophage VEH vs TAD") +
  theme(plot.title = element_text(size = 20, face = "bold"))

vol_plot

a=data.frame(Tcell_4v5DE)
if(a$padj<0.05){
    a$new_category="nontrivial"
} else {
    a$new_category="ns"
}
get_point_color <- function(x) {
  dplyr::case_when(x <= 0.05 ~ "black",
                   x > 0.05 ~ "red")
}
options(repr.plot.width =15, repr.plot.height =15)

vol_plot <- 
  ggplot(data=a,aes(x = log2FoldChange,
             y = -log10(padj),colour = get_point_color(a$padj))) + 
  geom_point(size=2) +
  geom_text(size=6,aes(label=ifelse(padj<0.05,as.character(rownames(a)),'')),hjust=0,vjust=0,) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  NoLegend() +
  ggtitle("vol_plot for Tcell MLIST vs MLIST.TAD") +
  theme(plot.title = element_text(size = 20, face = "bold"))

vol_plot

#### Heatmap for Granulocyte
## 1vs2
dfAnnot <- ONEvTWOmeta_Gra
HA <- HeatmapAnnotation(df = dfAnnot,
                        col=list(group_id =
                                   c('VEH'='blue',
                                     'TAD'='red')))


Heatmap(t(apply(datCombat[Gra_ONEvTWO_DEGenes, rownames(dfAnnot)],1,scale)),top_annotation = HA,
        clustering_distance_rows = 'pearson', 
        clustering_distance_columns = 'pearson',
        row_names_gp = gpar(fontsize = 6), column_title="Granulocyte 1v2 DE Genes")

### Heatmap
##AvB
dfAnnot <- AvBmeta[,4:6]
HA <- HeatmapAnnotation(df = dfAnnot,
                        col=list(Arm =
                                   c('A'='blue',
                                     'B'='red')))


Heatmap(t(apply(datCombat[AvBDEGenes, rownames(dfAnnot)],1,scale)),top_annotation = HA,
        clustering_distance_rows = 'pearson', 
        clustering_distance_columns = 'pearson',
        row_names_gp = gpar(fontsize = 6), column_title="CD8 AvB DE Genes")

##BvC
# dfAnnot <- BvCmeta[,4:6]
# HA <- HeatmapAnnotation(df = dfAnnot,
#                         col=list(Arm =
#                                    c('B'='blue',
#                                      'C'='red')))
# 
# 
# Heatmap(t(apply(datCombat[BvCDEGenes, rownames(dfAnnot)],1,scale)),top_annotation = HA,
#         clustering_distance_rows = 'pearson', 
#         clustering_distance_columns = 'pearson',
#         row_names_gp = gpar(fontsize = 6), column_title="CD8 BvC DE Genes")

##AvC
dfAnnot <- AvCmeta[,4:6]
HA <- HeatmapAnnotation(df = dfAnnot,
                        col=list(Arm =
                                   c('A'='blue',
                                     'C'='red')))


Heatmap(t(apply(datCombat[AvCDEGenes, rownames(dfAnnot)],1,scale)),top_annotation = HA,
        clustering_distance_rows = 'pearson', 
        clustering_distance_columns = 'pearson',
        row_names_gp = gpar(fontsize = 4), column_title="CD8 AvC DE Genes")


##All DE genes plotted on all patients -- genes grouped by which comparison they came from
dfAnnot <- meta[,4:6]
HA <- HeatmapAnnotation(df = dfAnnot,
                        col=list(Arm =
                                   c('A'='blue',
                                     'B' = 'yellow',
                                     'C'='red'),
                                 OS=
                                   c('<2years'='lightgreen',
                                     '>2years' = 'lightblue')))

AvBnAvCgenes <- intersect(AvBDEGenes,AvCDEGenes)
AvBDEGenes_filt <- AvBDEGenes[!AvBDEGenes %in% AvBnAvCgenes]
AvCDEGenes_filt <- AvCDEGenes[!AvCDEGenes %in% AvBnAvCgenes]

AvB <- cbind(AvBDEGenes_filt, rep("AvB", length(AvBDEGenes_filt)))
BvC <- cbind(BvCDEGenes, rep("BvC", length(BvCDEGenes)))
AvC <- cbind(AvCDEGenes_filt, rep("AvC", length(AvCDEGenes_filt)))
AvBnAvC <- cbind(AvBnAvCgenes, rep("AvB & AvC", length(AvBnAvCgenes)))

Allgenes <- as.data.frame(rbind(AvB,BvC,AvC,AvBnAvC))
colnames(Allgenes) <- c("Genes","Comparison")

Heatmap(t(apply(datCombat[Allgenes$Genes, rownames(dfAnnot)],1,scale)),top_annotation = HA,
        clustering_distance_rows = 'pearson', 
        clustering_distance_columns = 'pearson',
        row_names_gp = gpar(fontsize = 3), column_title="CD8 All DE Genes on all patients",
        row_split = Allgenes$Comparison, border = TRUE)






### HALLMARK statistics

## Print results in table and save
fwrite(as.data.frame(Gra_DE$ONEvTWOgs$HALLMARK), file="Granulocyte_HALLMARK_1v2.csv")
fwrite(as.data.frame(Gra_DE$FOURvFIVEgs$HALLMARK), file="Granulocyte_HALLMARK_4v5.csv")
fwrite(as.data.frame(Mac_DE$ONEvTWOgs$HALLMARK), file="Macrophage_HALLMARK_1v2.csv")
fwrite(as.data.frame(Mac_DE$FOURvFIVEgs$HALLMARK), file="Macrophage_HALLMARK_4v5.csv")
fwrite(as.data.frame(Tcell_DE$ONEvTWOgs$HALLMARK), file="Tcell_HALLMARK_1v2.csv")
fwrite(as.data.frame(Tcell_DE$FOURvFIVEgs$HALLMARK), file="Tcell_HALLMARK_4v5.csv")

#datatable(Gra_DE$ONEvTWOgs$HALLMARK)

#datatable(as.matrix(Gra_DE$one_twoRES))

Gra_DE$FOURvFIVEgs$HALLMARK$adjPvalue <- ifelse(Gra_DE$FOURvFIVEgs$HALLMARK$padj <= 0.1, "significant", "non-significant")
Gra_Sig_4v5ResHALLMARK <- Gra_DE$FOURvFIVEgs$HALLMARK[Gra_DE$FOURvFIVEgs$HALLMARK$adjPvalue=='significant']
Gra_Sig_4v5ResHALLMARK$adjPvalue <- ifelse(Gra_Sig_4v5ResHALLMARK$padj <= 0.05, "significant", "non-significant")

if (nrow(Gra_Sig_4v5ResHALLMARK) == 0) {
  print("No Significant pathways with pAdj < 0.05")
} else {
  ggplot(Gra_Sig_4v5ResHALLMARK, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
    geom_bar(width = 0.3, stat = "identity") +
    scale_fill_manual(values = c("significant" = "red", "non-significant" = "grey")) +
    theme(axis.text.y = element_text(size=5),
          plot.title = element_text(size=8, face="bold"),
          plot.margin = unit(c(1, 1, 1, 1), "cm"),
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white")) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="Significant (pAdj<0.05) Hallmark pathways EA from 4v5 in Granulocyte") + 
    theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

## Granulocyte
# 1v2
Gra_DE$ONEvTWOgs$HALLMARK$adjPvalue <- ifelse(Gra_DE$ONEvTWOgs$HALLMARK$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Gra_Sig_1v2ResHALLMARK <- Gra_DE$ONEvTWOgs$HALLMARK[padj < 0.05]
if (nrow(Gra_Sig_1v2ResHALLMARK)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Gra_Sig_1v2ResHALLMARK, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_bar(width = 0.3, stat = "identity") +  # Set the width parameter here
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5),
        plot.title = element_text(size=8, face="bold"),  # Adjust title size
        plot.margin = unit(c(1, 1, 1, 1), "cm"),       # Adjust margins (top, right, bottom, left)
        plot.background = element_rect(fill = "white"), # Set plot background color
        panel.background = element_rect(fill = "white")) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) Hallmark pathways EA from 1v2 in Granulocyte") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Gra_DE$FOURvFIVEgs$HALLMARK$adjPvalue <- ifelse(Gra_DE$FOURvFIVEgs$HALLMARK$padj <= 0.05, "significant", "non-significant")
#cols <- c("non-significant" = "grey", "significant" = "red")

Gra_Sig_4v5ResHALLMARK <- Gra_DE$FOURvFIVEgs$HALLMARK[padj < 0.05]
if (nrow(Gra_Sig_4v5ResHALLMARK)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Gra_Sig_4v5ResHALLMARK, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_bar(width = 0.3, stat = "identity") +  # Set the width parameter here
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5),
        plot.title = element_text(size=8, face="bold"),  # Adjust title size
        plot.margin = unit(c(1, 1, 1, 1), "cm"),       # Adjust margins (top, right, bottom, left)
        plot.background = element_rect(fill = "white"), # Set plot background color
        panel.background = element_rect(fill = "white")) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) Hallmark pathways EA from 4v5 in Granulocyte") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

## Macrophage
# 1v2
Mac_DE$ONEvTWOgs$HALLMARK$adjPvalue <- ifelse(Mac_DE$ONEvTWOgs$HALLMARK$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Mac_Sig_1v2ResHALLMARK <- Mac_DE$ONEvTWOgs$HALLMARK[padj < 0.05]
if (nrow(Mac_Sig_1v2ResHALLMARK)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Mac_Sig_1v2ResHALLMARK, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_bar(width = 0.7, stat = "identity") +  # Set the width parameter here
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) Hallmark pathways EA from 1v2 in Macrophage") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Mac_DE$FOURvFIVEgs$HALLMARK$adjPvalue <- ifelse(Mac_DE$FOURvFIVEgs$HALLMARK$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Mac_Sig_4v5ResHALLMARK <- Mac_DE$FOURvFIVEgs$HALLMARK[padj < 0.05]
if (nrow(Mac_Sig_4v5ResHALLMARK)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Mac_Sig_4v5ResHALLMARK, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_bar(width = 0.7, stat = "identity") +  # Set the width parameter here
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) Hallmark pathways EA from 4v5 in Macrophage") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

backup <- Tcell_DE$ONEvTWOgs$HALLMARK[Tcell_DE$ONEvTWOgs$HALLMARK$pathway=="HALLMARK_IL2_STAT5_SIGNALING",]
backup2 <- Tcell_DE$ONEvTWOgs$HALLMARK[Tcell_DE$ONEvTWOgs$HALLMARK$pathway=="HALLMARK_INFLAMMATORY_RESPONSE",]

ggplot(backup2, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) Hallmark pathways EA from 1v2 in Macrophage") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))

## Tcell
# 1v2
Tcell_DE$ONEvTWOgs$HALLMARK$adjPvalue <- ifelse(Tcell_DE$ONEvTWOgs$HALLMARK$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Tcell_Sig_1v2ResHALLMARK <- Tcell_DE$ONEvTWOgs$HALLMARK[padj < 0.05]
if (nrow(Tcell_Sig_1v2ResHALLMARK)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Tcell_Sig_1v2ResHALLMARK, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) Hallmark pathways EA from 1v2 in Macrophage") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Tcell_DE$FOURvFIVEgs$HALLMARK$adjPvalue <- ifelse(Tcell_DE$FOURvFIVEgs$HALLMARK$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Tcell_Sig_4v5ResHALLMARK <- Tcell_DE$FOURvFIVEgs$HALLMARK[padj < 0.05]
if (nrow(Tcell_Sig_4v5ResHALLMARK)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Tcell_Sig_4v5ResHALLMARK, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) Hallmark pathways EA from 4v5 in Tcell") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

### KEGG statistics
# 1v2
## Print results in table and save
fwrite(as.data.frame(Gra_DE$ONEvTWOgs$KEGG), file="Granulocyte_KEGG_1v2Treatment.csv")
fwrite(as.data.frame(Gra_DE$FOURvFIVEgs$KEGG), file="Granulocyte_KEGG_4v5Treatment.csv")
fwrite(as.data.frame(Mac_DE$ONEvTWOgs$KEGG), file="Macrophage_KEGG_1v2Treatment.csv")
fwrite(as.data.frame(Mac_DE$FOURvFIVEgs$KEGG), file="Macrophage_KEGG_4v5Treatment.csv")
fwrite(as.data.frame(Tcell_DE$ONEvTWOgs$KEGG), file="Tcell_KEGG_1v2Treatment.csv")
fwrite(as.data.frame(Tcell_DE$FOURvFIVEgs$KEGG), file="Tcell_KEGG_4v5Treatment.csv")

## Granulocyte
#1v2
Gra_DE$ONEvTWOgs$KEGG$adjPvalue <- ifelse(Gra_DE$ONEvTWOgs$KEGG$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Gra_Sig_1v2ResKEGG <- Gra_DE$ONEvTWOgs$KEGG[padj < 0.05]
if (nrow(Gra_Sig_1v2ResKEGG)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Gra_Sig_1v2ResKEGG, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) KEGG pathways Enrichment Score from Treatment (1v2) in Granulocyte Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Gra_DE$FOURvFIVEgs$KEGG$adjPvalue <- ifelse(Gra_DE$FOURvFIVEgs$KEGG$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Gra_Sig_4v5ResKEGG <- Gra_DE$FOURvFIVEgs$KEGG[padj < 0.05]
if (nrow(Gra_Sig_4v5ResKEGG)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Gra_Sig_4v5ResKEGG, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) KEGG pathways Enrichment Score from Treatment (4v5) in Granulocyte Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

## Macrophage
#1v2
Mac_DE$ONEvTWOgs$KEGG$adjPvalue <- ifelse(Mac_DE$ONEvTWOgs$KEGG$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Mac_Sig_1v2ResKEGG <- Mac_DE$ONEvTWOgs$KEGG[padj < 0.05]
if (nrow(Mac_Sig_1v2ResKEGG)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Mac_Sig_1v2ResKEGG, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) KEGG pathways Enrichment Score from Treatment (1v2) in Macrophage Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}
#4v5
Mac_DE$FOURvFIVEgs$KEGG$adjPvalue <- ifelse(Mac_DE$FOURvFIVEgs$KEGG$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Mac_Sig_4v5ResKEGG <- Mac_DE$FOURvFIVEgs$KEGG[padj < 0.05]
if (nrow(Mac_Sig_4v5ResKEGG)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Mac_Sig_4v5ResKEGG, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) KEGG pathways Enrichment Score from Treatment (4v5) in Macrophage Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

## Tcell
#1v2

options(repr.plot.width =20, repr.plot.height =3)

Tcell_DE$ONEvTWOgs$KEGG$adjPvalue <- ifelse(Tcell_DE$ONEvTWOgs$KEGG$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Tcell_Sig_1v2ResKEGG <- Tcell_DE$ONEvTWOgs$KEGG[padj < 0.05]
if (nrow(Tcell_Sig_1v2ResKEGG)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Tcell_Sig_1v2ResKEGG, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) KEGG pathways Enrichment Score from Treatment (1v2) in Tcell Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Tcell_DE$FOURvFIVEgs$KEGG$adjPvalue <- ifelse(Tcell_DE$FOURvFIVEgs$KEGG$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Tcell_Sig_4v5ResKEGG <- Tcell_DE$FOURvFIVEgs$KEGG[padj < 0.05]
if (nrow(Tcell_Sig_4v5ResKEGG)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Tcell_Sig_4v5ResKEGG, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) KEGG pathways Enrichment Score from Treatment (4v5) in Tcell Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

### REACTOME statistics
## Print results in table and save
fwrite(as.data.frame(Gra_DE$ONEvTWOgs$REACTOME), file="Granulocyte_REACTOME_1v2Treatment.csv")
fwrite(as.data.frame(Gra_DE$FOURvFIVEgs$REACTOME), file="Granulocyte_REACTOME_4v5Treatment.csv")
fwrite(as.data.frame(Mac_DE$ONEvTWOgs$REACTOME), file="Macrophage_REACTOME_1v2Treatment.csv")
fwrite(as.data.frame(Mac_DE$FOURvFIVEgs$REACTOME), file="Macrophage_REACTOME_4v5Treatment.csv")
fwrite(as.data.frame(Tcell_DE$ONEvTWOgs$REACTOME), file="Tcell_REACTOME_1v2Treatment.csv")
fwrite(as.data.frame(Tcell_DE$FOURvFIVEgs$REACTOME), file="Tcell_REACTOME_4v5Treatment.csv")

## Macrophage
# 1v2
options(repr.plot.width =20, repr.plot.height =9)

Mac_DE$ONEvTWOgs$REACTOME$adjPvalue <- ifelse(Mac_DE$ONEvTWOgs$REACTOME$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Mac_Sig_1v2ResREACTOME <- Mac_DE$ONEvTWO$REACTOME[padj < 0.05]
if (nrow(Mac_Sig_1v2ResREACTOME)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Mac_Sig_1v2ResREACTOME, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) REACTOME pathways Enrichment Score from Treatment (4v5) in Macrophage Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Mac_DE$FOURvFIVEgs$REACTOME$adjPvalue <- ifelse(Mac_DE$FOURvFIVEgs$REACTOME$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Mac_Sig_4v5ResREACTOME <- Mac_DE$FOURvFIVEgs$REACTOME[padj < 0.05]
if (nrow(Mac_Sig_4v5ResREACTOME)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Mac_Sig_4v5ResREACTOME, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) REACTOME pathways Enrichment Score from Treatment (4v5) in Macrophage Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

## Granulocyte
#1v2
Gra_DE$ONEvTWOgs$REACTOME$adjPvalue <- ifelse(Gra_DE$ONEvTWOgs$REACTOME$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Gra_Sig_1v2ResREACTOME <- Gra_DE$ONEvTWOgs$REACTOME[padj < 0.05]
if (nrow(Gra_Sig_1v2ResREACTOME)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Gra_Sig_1v2ResREACTOME, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) REACTOME pathways Enrichment Score from Treatment (1v2) in Tcell Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Gra_DE$FOURvFIVEgs$REACTOME$adjPvalue <- ifelse(Gra_DE$FOURvFIVEgs$REACTOME$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Gra_Sig_4v5ResREACTOME <- Gra_DE$FOURvFIVEgs$REACTOME[padj < 0.05]
if (nrow(Gra_Sig_4v5ResREACTOME)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Gra_Sig_4v5ResREACTOME, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) REACTOME pathways Enrichment Score from Treatment (4v5) in Granulocyte Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

## Tcell
#1v2
Tcell_DE$ONEvTWOgs$REACTOME$adjPvalue <- ifelse(Tcell_DE$ONEvTWOgs$REACTOME$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Tcell_Sig_1v2ResREACTOME <- Tcell_DE$ONEvTWOgs$REACTOME[padj < 0.05]
if (nrow(Tcell_Sig_1v2ResREACTOME)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Tcell_Sig_1v2ResREACTOME, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) REACTOME pathways Enrichment Score from Treatment (1v2) in Tcell Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Tcell_DE$FOURvFIVEgs$REACTOME$adjPvalue <- ifelse(Tcell_DE$FOURvFIVEgs$REACTOME$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Tcell_Sig_4v5ResREACTOME <- Tcell_DE$FOURvFIVEgs$REACTOME[padj < 0.05]
if (nrow(Tcell_Sig_4v5ResREACTOME)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Tcell_Sig_4v5ResREACTOME, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) REACTOME pathways Enrichment Score from Treatment (4v5) in Tcell Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

### IMMUNESIG statistics
## Print results in table and save
fwrite(as.data.frame(Gra_DE$ONEvTWOgs$IMMUNESIG), file="Granulocyte_IMMUNESIG_1v2Treatment.csv")
fwrite(as.data.frame(Gra_DE$FOURvFIVEgs$IMMUNESIG), file="Granulocyte_IMMUNESIG_4v5Treatment.csv")
fwrite(as.data.frame(Mac_DE$ONEvTWOgs$IMMUNESIG), file="Macrophage_IMMUNESIG_1v2Treatment.csv")
fwrite(as.data.frame(Mac_DE$FOURvFIVEgs$IMMUNESIG), file="Macrophage_IMMUNESIG_4v5Treatment.csv")
fwrite(as.data.frame(Tcell_DE$ONEvTWOgs$IMMUNESIG), file="Tcell_IMMUNESIG_1v2Treatment.csv")
fwrite(as.data.frame(Tcell_DE$FOURvFIVEgs$IMMUNESIG), file="Tcell_IMMUNESIG_4v5Treatment.csv")

## Granulocyte
#1v2
options(repr.plot.width =25, repr.plot.height =100)
Gra_DE$ONEvTWOgs$IMMUNESIG$adjPvalue <- ifelse(Gra_DE$ONEvTWOgs$IMMUNESIG$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Gra_Sig_1v2ResIMMUNESIG <- Gra_DE$ONEvTWOgs$IMMUNESIG[padj < 0.05]
if (nrow(Gra_Sig_1v2ResIMMUNESIG)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Gra_Sig_1v2ResIMMUNESIG, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) IMMUNESIG pathways Enrichment Score from Treatment (1v2) in Granulocyte Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Gra_DE$FOURvFIVEgs$IMMUNESIG$adjPvalue <- ifelse(Gra_DE$FOURvFIVEgs$IMMUNESIG$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Gra_Sig_4v5ResIMMUNESIG <- Gra_DE$FOURvFIVEgs$IMMUNESIG[padj < 0.05]
if (nrow(Gra_Sig_4v5ResIMMUNESIG)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Gra_Sig_4v5ResIMMUNESIG, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) IMMUNESIG pathways Enrichment Score from Treatment (4v5) in Granulocyte Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

## Macrophage
#1v2
options(repr.plot.width =25, repr.plot.height =20)

Mac_DE$ONEvTWOgs$IMMUNESIG$adjPvalue <- ifelse(Mac_DE$ONEvTWOgs$IMMUNESIG$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Mac_Sig_1v2ResIMMUNESIG <- Mac_DE$ONEvTWOgs$IMMUNESIG[padj < 0.05]
if (nrow(Mac_Sig_1v2ResIMMUNESIG)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Mac_Sig_1v2ResIMMUNESIG, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) IMMUNESIG pathways Enrichment Score from Treatment (1v2) in Macrophage Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Mac_DE$FOURvFIVEgs$IMMUNESIG$adjPvalue <- ifelse(Mac_DE$FOURvFIVEgs$IMMUNESIG$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Mac_Sig_4v5ResIMMUNESIG <- Mac_DE$FOURvFIVEgs$IMMUNESIG[padj < 0.05]
if (nrow(Mac_Sig_4v5ResIMMUNESIG)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Mac_Sig_4v5ResIMMUNESIG, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) IMMUNESIG pathways Enrichment Score from Treatment (4v5) in Macrophage Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

## Tcell
#1v2
options(repr.plot.width =25, repr.plot.height =20)

Tcell_DE$ONEvTWOgs$IMMUNESIG$adjPvalue <- ifelse(Tcell_DE$ONEvTWOgs$IMMUNESIG$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Tcell_Sig_1v2ResIMMUNESIG <- Tcell_DE$ONEvTWOgs$IMMUNESIG[padj < 0.05]
if (nrow(Tcell_Sig_1v2ResIMMUNESIG)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Tcell_Sig_1v2ResIMMUNESIG, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) IMMUNESIG pathways Enrichment Score from Treatment (1v2) in Macrophage Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Tcell_DE$FOURvFIVEgs$IMMUNESIG$adjPvalue <- ifelse(Tcell_DE$FOURvFIVEgs$IMMUNESIG$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Tcell_Sig_4v5ResIMMUNESIG <- Tcell_DE$FOURvFIVEgs$IMMUNESIG[padj < 0.05]
if (nrow(Tcell_Sig_4v5ResIMMUNESIG)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Tcell_Sig_4v5ResIMMUNESIG, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) IMMUNESIG pathways Enrichment Score from Treatment (4v5) in Macrophage Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

### GOBP statistics
## Print results in table and save
fwrite(as.data.frame(Gra_DE$ONEvTWOgs$GOBP), file="Granulocyte_GOBP_1v2Treatment.csv")
fwrite(as.data.frame(Gra_DE$FOURvFIVEgs$GOBP), file="Granulocyte_GOBP_4v5Treatment.csv")
fwrite(as.data.frame(Mac_DE$ONEvTWOgs$GOBP), file="Macrophage_GOBP_1v2Treatment.csv")
fwrite(as.data.frame(Mac_DE$FOURvFIVEgs$GOBP), file="Macrophage_GOBP_4v5Treatment.csv")
fwrite(as.data.frame(Tcell_DE$ONEvTWOgs$GOBP), file="Tcell_GOBP_1v2Treatment.csv")
fwrite(as.data.frame(Tcell_DE$FOURvFIVEgs$GOBP), file="Tcell_GOBP_4v5Treatment.csv")

## Granulocyte
#1v2
options(repr.plot.width =20, repr.plot.height =20)

Gra_DE$ONEvTWOgs$GOBP$adjPvalue <- ifelse(Gra_DE$ONEvTWOgs$GOBP$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Gra_Sig_1v2ResGOBP <- Gra_DE$ONEvTWOgs$GOBP[padj < 0.05]
if (nrow(Gra_Sig_1v2ResGOBP)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Gra_Sig_1v2ResGOBP, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) GOBP pathways Enrichment Score from Treatment (1v2) in Granulocyte Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Gra_DE$FOURvFIVEgs$GOBP$adjPvalue <- ifelse(Gra_DE$FOURvFIVEgs$GOBP$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Gra_Sig_4v5ResGOBP <- Gra_DE$FOURvFIVEgs$GOBP[padj < 0.05]
if (nrow(Gra_Sig_4v5ResGOBP)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Gra_Sig_4v5ResGOBP, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) GOBP pathways Enrichment Score from Treatment (4v5) in Granulocyte Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

## Macrophage
#1v2
options(repr.plot.width =20, repr.plot.height =7)

Mac_DE$ONEvTWOgs$GOBP$adjPvalue <- ifelse(Mac_DE$ONEvTWOgs$GOBP$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Mac_Sig_1v2ResGOBP <- Mac_DE$ONEvTWOgs$GOBP[padj < 0.05]
if (nrow(Mac_Sig_1v2ResGOBP)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Mac_Sig_1v2ResGOBP, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) GOBP pathways Enrichment Score from Treatment (4v5) in Macrophage Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Mac_DE$FOURvFIVEgs$GOBP$adjPvalue <- ifelse(Mac_DE$FOURvFIVEgs$GOBP$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Mac_Sig_4v5ResGOBP <- Mac_DE$FOURvFIVEgs$GOBP[padj < 0.05]
if (nrow(Mac_Sig_4v5ResGOBP)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Mac_Sig_4v5ResGOBP, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) GOBP pathways Enrichment Score from Treatment (4v5) in Macrophage Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

## Tcell
#1v2
options(repr.plot.width =20, repr.plot.height =5)

Tcell_DE$ONEvTWOgs$GOBP$adjPvalue <- ifelse(Tcell_DE$ONEvTWOgs$GOBP$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Tcell_Sig_1v2ResGOBP <- Tcell_DE$ONEvTWOgs$GOBP[padj < 0.05]
if (nrow(Tcell_Sig_1v2ResGOBP)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Tcell_Sig_1v2ResGOBP, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) GOBP pathways Enrichment Score from Treatment (1v2) in Tcell") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Tcell_DE$FOURvFIVEgs$GOBP$adjPvalue <- ifelse(Tcell_DE$FOURvFIVEgs$GOBP$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Tcell_Sig_4v5ResGOBP <- Tcell_DE$FOURvFIVEgs$GOBP[padj < 0.05]
if (nrow(Tcell_Sig_4v5ResGOBP)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Tcell_Sig_4v5ResGOBP, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) GOBP pathways Enrichment Score from Treatment (4v5) in ") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

### GOMF statistics
## Print results in table and save
fwrite(as.data.frame(Gra_DE$ONEvTWOgs$GOMF), file="Granulocyte_GOMF_1v2Treatment.csv")
fwrite(as.data.frame(Gra_DE$FOURvFIVEgs$GOMF), file="Granulocyte_GOMF_4v5Treatment.csv")
fwrite(as.data.frame(Mac_DE$ONEvTWOgs$GOMF), file="Macrophage_GOMF_1v2Treatment.csv")
fwrite(as.data.frame(Mac_DE$FOURvFIVEgs$GOMF), file="Macrophage_GOMF_4v5Treatment.csv")
fwrite(as.data.frame(Tcell_DE$ONEvTWOgs$GOMF), file="Tcell_GOMF_1v2Treatment.csv")
fwrite(as.data.frame(Tcell_DE$FOURvFIVEgs$GOMF), file="Tcell_GOMF_4v5Treatment.csv")

## Granulocyte
#1v2
options(repr.plot.width =20, repr.plot.height =20)

Gra_DE$ONEvTWOgs$GOMF$adjPvalue <- ifelse(Gra_DE$ONEvTWOgs$GOMF$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Gra_Sig_1v2ResGOMF <- Gra_DE$ONEvTWOgs$GOMF[padj < 0.05]
if (nrow(Gra_Sig_1v2ResGOMF)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Gra_Sig_1v2ResGOMF, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) GOMF pathways Enrichment Score from Treatment (1v2) in Granulocyte Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Gra_DE$FOURvFIVEgs$GOMF$adjPvalue <- ifelse(Gra_DE$FOURvFIVEgs$GOMF$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Gra_Sig_4v5ResGOMF <- Gra_DE$FOURvFIVEgs$GOMF[padj < 0.05]
if (nrow(Gra_Sig_4v5ResGOMF)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Gra_Sig_4v5ResGOMF, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) GOMF pathways Enrichment Score from Treatment (4v5) in Granulocyte Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

## Macrophage
#1v2
options(repr.plot.width =10, repr.plot.height =2)

Mac_DE$ONEvTWOgs$GOMF$adjPvalue <- ifelse(Mac_DE$ONEvTWOgs$GOMF$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Mac_Sig_1v2ResGOMF <- Mac_DE$ONEvTWOgs$GOMF[padj < 0.05]
if (nrow(Mac_Sig_1v2ResGOMF)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Mac_Sig_1v2ResGOMF, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) GOMF pathways Enrichment Score from Treatment (1v2) in Macrophage Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Mac_DE$FOURvFIVEgs$GOMF$adjPvalue <- ifelse(Mac_DE$FOURvFIVEgs$GOMF$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Mac_Sig_4v5ResGOMF <- Mac_DE$FOURvFIVEgs$GOMF[padj < 0.05]
if (nrow(Mac_Sig_4v5ResGOMF)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Mac_Sig_4v5ResGOMF, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) GOMF pathways Enrichment Score from Treatment (4v5) in Macrophage Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

## Tcell
#1v2
options(repr.plot.width =10, repr.plot.height =2)

Tcell_DE$ONEvTWOgs$GOMF$adjPvalue <- ifelse(Tcell_DE$ONEvTWOgs$GOMF$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Tcell_Sig_1v2ResGOMF <- Tcell_DE$ONEvTWOgs$GOMF[padj < 0.05]
if (nrow(Tcell_Sig_1v2ResGOMF)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Tcell_Sig_1v2ResGOMF, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) GOMF pathways Enrichment Score from Treatment (1v2) in Macrophage Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Tcell_DE$FOURvFIVEgs$GOMF$adjPvalue <- ifelse(Tcell_DE$FOURvFIVEgs$GOMF$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Tcell_Sig_4v5ResGOMF <- Tcell_DE$FOURvFIVEgs$GOMF[padj < 0.05]
if (nrow(Tcell_Sig_4v5ResGOMF)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Tcell_Sig_4v5ResGOMF, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) GOMF pathways Enrichment Score from Treatment (4v5) in Macrophage Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

### NABA statistics
## Print results in table and save
fwrite(as.data.frame(Gra_DE$ONEvTWOgs$NABA), file="Granulocyte_NABA_1v2Treatment.csv")
fwrite(as.data.frame(Gra_DE$FOURvFIVEgs$NABA), file="Granulocyte_NABA_4v5Treatment.csv")
fwrite(as.data.frame(Mac_DE$ONEvTWOgs$NABA), file="Macrophage_NABA_1v2Treatment.csv")
fwrite(as.data.frame(Mac_DE$FOURvFIVEgs$NABA), file="Macrophage_NABA_4v5Treatment.csv")
fwrite(as.data.frame(Tcell_DE$ONEvTWOgs$NABA), file="Tcell_NABA_1v2Treatment.csv")
fwrite(as.data.frame(Tcell_DE$FOURvFIVEgs$NABA), file="Tcell_NABA_4v5Treatment.csv")

## Granulocyte
#1v2
options(repr.plot.width =20, repr.plot.height =7)

Gra_DE$ONEvTWOgs$NABA$adjPvalue <- ifelse(Gra_DE$ONEvTWOgs$NABA$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Gra_Sig_1v2ResNABA <- Gra_DE$ONEvTWOgs$NABA[padj < 0.05]
if (nrow(Gra_Sig_1v2ResNABA)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Gra_Sig_1v2ResNABA, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) NABA pathways Enrichment Score from Treatment (1v2) in Granulocyte Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Gra_DE$FOURvFIVEgs$NABA$adjPvalue <- ifelse(Gra_DE$FOURvFIVEgs$NABA$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Gra_Sig_4v5ResNABA <- Gra_DE$FOURvFIVEgs$NABA[padj < 0.05]
if (nrow(Gra_Sig_4v5ResNABA)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Gra_Sig_4v5ResNABA, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) NABA pathways Enrichment Score from Treatment (4v5) in Granulocyte Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

## Macrophage
#1v2
options(repr.plot.width =20, repr.plot.height =7)

Mac_DE$ONEvTWOgs$NABA$adjPvalue <- ifelse(Mac_DE$ONEvTWOgs$NABA$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Mac_Sig_1v2ResNABA <- Mac_DE$ONEvTWOgs$NABA[padj < 0.05]
if (nrow(Mac_Sig_1v2ResNABA)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Mac_Sig_1v2ResNABA, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) NABA pathways Enrichment Score from Treatment (1v2) in Granulocyte Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Mac_DE$FOURvFIVEgs$NABA$adjPvalue <- ifelse(Mac_DE$FOURvFIVEgs$NABA$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Mac_Sig_4v5ResNABA <- Mac_DE$FOURvFIVEgs$NABA[padj < 0.05]
if (nrow(Mac_Sig_4v5ResNABA)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Mac_Sig_4v5ResNABA, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) NABA pathways Enrichment Score from Treatment (4v5) in Granulocyte Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

## Tcell
#1v2
options(repr.plot.width =20, repr.plot.height =7)

Tcell_DE$ONEvTWOgs$NABA$adjPvalue <- ifelse(Tcell_DE$ONEvTWOgs$NABA$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Tcell_Sig_1v2ResNABA <- Tcell_DE$ONEvTWOgs$NABA[padj < 0.05]
if (nrow(Tcell_Sig_1v2ResNABA)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Tcell_Sig_1v2ResNABA, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) NABA pathways Enrichment Score from Treatment (1v2) in T Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}

#4v5
Tcell_DE$FOURvFIVEgs$NABA$adjPvalue <- ifelse(Tcell_DE$FOURvFIVEgs$NABA$padj  <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

Tcell_Sig_4v5ResNABA <- Tcell_DE$FOURvFIVEgs$NABA[padj < 0.05]
if (nrow(Tcell_Sig_4v5ResNABA)==0){ print("No Significant pathways with pAdj < 0.05"
)} else {
  ggplot(Tcell_Sig_4v5ResNABA, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.y = element_text(size=5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
  title="Significant (pAdj<0.05) NABA pathways Enrichment Score from Treatment (4v5) in T Cells") + theme(plot.title.position = "plot", plot.title = element_text(size=8, face="bold"))
}



### extract top genes from important pathway sets

##### 1v2 Granulocyte

Gra1v2 <- data.frame(baseMean=NULL,log2FoldChange=NULL,lfcSE=NULL,stat=NULL,pvalue=NULL,padj=NULL,pathway=NULL,celltype=NULL,comparison=NULL,mean_VEH=NULL,mean_TAD=NULL)

### HALLMARK Granulocyte 1v2
##  HALLMARK_INFLAMMATORY_RESPONSE
top_gene1=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$HALLMARK))
top_gene1=top_gene1[top_gene1$pathway=='HALLMARK_INFLAMMATORY_RESPONSE',]
top_gene1=unlist(top_gene1[,8])
HALL_INF_RES=as.data.frame(Gra_1v2DE)
HALL_INF_RES=HALL_INF_RES[(row.names(HALL_INF_RES) %in% top_gene1),]
HALL_INF_RES=HALL_INF_RES[HALL_INF_RES$pvalue<0.05,]
HALL_INF_RES$pathway="HALLMARK_INFLAMMATORY_RESPONSE"
HALL_INF_RES$celltype="Granulocyte"
HALL_INF_RES$comparison="VEH and TAD"

### HALLMARK Granulocyte 1v2
##  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
top_gene2=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$HALLMARK))
top_gene2=top_gene2[top_gene2$pathway=='HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',]
top_gene2=unlist(top_gene2[,8])
HALL_EPI_MES_TRANS=as.data.frame(Gra_1v2DE)
HALL_EPI_MES_TRANS=HALL_EPI_MES_TRANS[(row.names(HALL_EPI_MES_TRANS) %in% top_gene2),]
HALL_EPI_MES_TRANS=HALL_EPI_MES_TRANS[HALL_EPI_MES_TRANS$pvalue<0.05,]
HALL_EPI_MES_TRANS$pathway="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
HALL_EPI_MES_TRANS$celltype="Granulocyte"
HALL_EPI_MES_TRANS$comparison="VEH and TAD"

### HALLMARK Granulocyte 1v2
##  HALLMARK_IL2_STAT5_SIGNALING
top_gene3=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$HALLMARK))
top_gene3=top_gene3[top_gene3$pathway=='HALLMARK_IL2_STAT5_SIGNALING',]
top_gene3=unlist(top_gene3[,8])
HALL_IL2=as.data.frame(Gra_1v2DE)
HALL_IL2=HALL_IL2[(row.names(HALL_IL2) %in% top_gene3),]
HALL_IL2=HALL_IL2[HALL_IL2$pvalue<0.05,]
HALL_IL2$pathway="HALLMARK_IL2_STAT5_SIGNALING"
HALL_IL2$celltype="Granulocyte"
HALL_IL2$comparison="VEH and TAD"

### KEGG Granulocyte 1v2
##  KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION
top_gene4=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$KEGG))
top_gene4=top_gene4[top_gene4$pathway=='KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION',]
top_gene4=unlist(top_gene4[,8])
KEGG_CYTO_REC=as.data.frame(Gra_1v2DE)
KEGG_CYTO_REC=KEGG_CYTO_REC[(row.names(KEGG_CYTO_REC) %in% top_gene4),]
KEGG_CYTO_REC=KEGG_CYTO_REC[,]
KEGG_CYTO_REC$pathway="KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION"
KEGG_CYTO_REC$celltype="Granulocyte"
KEGG_CYTO_REC$comparison="VEH and TAD"

### KEGG Granulocyte 1v2
##  KEGG_LYSOSOME
top_gene5=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$KEGG))
top_gene5=top_gene5[top_gene5$pathway=='KEGG_LYSOSOME',]
top_gene5=unlist(top_gene5[,8])
KEGG_LYSO=as.data.frame(Gra_1v2DE)
KEGG_LYSO=KEGG_LYSO[(row.names(KEGG_LYSO) %in% top_gene5),]
KEGG_LYSO=KEGG_LYSO[KEGG_LYSO$pvalue<0.05,]
KEGG_LYSO$pathway="KEGG_LYSOSOME"
KEGG_LYSO$celltype="Granulocyte"
KEGG_LYSO$comparison="VEH and TAD"

### REACTOME Granulocyte 1V2
##  REACTOME_INNATE_IMMUNE_SYSTEM
top_gene6=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$REACTOME))
top_gene6=top_gene6[top_gene6$pathway=='REACTOME_INNATE_IMMUNE_SYSTEM',]
top_gene6=unlist(top_gene6[,8])
REACTOME_INNATE=as.data.frame(Gra_1v2DE)
REACTOME_INNATE=REACTOME_INNATE[(row.names(REACTOME_INNATE) %in% top_gene6),]
REACTOME_INNATE=REACTOME_INNATE[REACTOME_INNATE$pvalue<0.05,]
REACTOME_INNATE$pathway="REACTOME_INNATE_IMMUNE_SYSTEM"
REACTOME_INNATE$celltype="Granulocyte"
REACTOME_INNATE$comparison="VEH and TAD"

### REACTOME Granulocyte 1V2
##  REACTOME_NEUTROPHIL_DEGRANULATION
top_gene7=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$REACTOME))
top_gene7=top_gene7[top_gene7$pathway=='REACTOME_NEUTROPHIL_DEGRANULATION',]
top_gene7=unlist(top_gene7[,8])
REACTOME_NEUTRO=as.data.frame(Gra_1v2DE)
REACTOME_NEUTRO=REACTOME_NEUTRO[(row.names(REACTOME_NEUTRO) %in% top_gene7),]
REACTOME_NEUTRO=REACTOME_NEUTRO[REACTOME_NEUTRO$pvalue<0.05,]
REACTOME_NEUTRO$pathway="REACTOME_NEUTROPHIL_DEGRANULATION"
REACTOME_NEUTRO$celltype="Granulocyte"
REACTOME_NEUTRO$comparison="VEH and TAD"

### REACTOME Granulocyte 1V2
##  REACTOME_COLLAGEN_FORMATION
top_gene9=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$REACTOME))
top_gene9=top_gene9[top_gene9$pathway=='REACTOME_COLLAGEN_FORMATION',]
top_gene9=unlist(top_gene9[,8])
REACTOME_COLLAGEN_FORM=as.data.frame(Gra_1v2DE)
REACTOME_COLLAGEN_FORM=REACTOME_COLLAGEN_FORM[(row.names(REACTOME_COLLAGEN_FORM) %in% top_gene9),]
REACTOME_COLLAGEN_FORM=REACTOME_COLLAGEN_FORM[REACTOME_COLLAGEN_FORM$pvalue<0.05,]
REACTOME_COLLAGEN_FORM$pathway="REACTOME_COLLAGEN_FORMATION"
REACTOME_COLLAGEN_FORM$celltype="Granulocyte"
REACTOME_COLLAGEN_FORM$comparison="VEH and TAD"

### REACTOME Granulocyte 1V2
##  REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES
top_gene10=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$REACTOME))
top_gene10=top_gene10[top_gene10$pathway=='REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES',]
top_gene10=unlist(top_gene10[,8])
REACTOME_ASSEMBLY_OTHER_STRUC=as.data.frame(Gra_1v2DE)
REACTOME_ASSEMBLY_OTHER_STRUC=REACTOME_ASSEMBLY_OTHER_STRUC[(row.names(REACTOME_ASSEMBLY_OTHER_STRUC) %in% top_gene10),]
REACTOME_ASSEMBLY_OTHER_STRUC=REACTOME_ASSEMBLY_OTHER_STRUC[REACTOME_ASSEMBLY_OTHER_STRUC$pvalue<0.05,]
REACTOME_ASSEMBLY_OTHER_STRUC$pathway="REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES"
REACTOME_ASSEMBLY_OTHER_STRUC$celltype="Granulocyte"
REACTOME_ASSEMBLY_OTHER_STRUC$comparison="VEH and TAD"


### REACTOME Granulocyte 1V2
##  REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION
top_gene12=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$REACTOME))
top_gene12=top_gene12[top_gene12$pathway=='REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION',]
top_gene12=unlist(top_gene12[,8])
REACTOME_EXTRA=as.data.frame(Gra_1v2DE)
REACTOME_EXTRA=REACTOME_EXTRA[(row.names(REACTOME_EXTRA) %in% top_gene12),]
REACTOME_EXTRA=REACTOME_EXTRA[REACTOME_EXTRA$pvalue<0.05,]
REACTOME_EXTRA$pathway="REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"
REACTOME_EXTRA$celltype="Granulocyte"
REACTOME_EXTRA$comparison="VEH and TAD"

### REACTOME Granulocyte 1V2
##  REACTOME_L1CAM_INTERACTIONS
top_gene13=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$REACTOME))
top_gene13=top_gene13[top_gene13$pathway=='REACTOME_L1CAM_INTERACTIONS',]
top_gene13=unlist(top_gene13[,8])
REACTOME_L1CAM=as.data.frame(Gra_1v2DE)
REACTOME_L1CAM=REACTOME_L1CAM[(row.names(REACTOME_L1CAM) %in% top_gene13),]
REACTOME_L1CAM=REACTOME_L1CAM[REACTOME_L1CAM$pvalue<0.05,]
REACTOME_L1CAM$pathway="REACTOME_L1CAM_INTERACTIONS"
REACTOME_L1CAM$celltype="Granulocyte"
REACTOME_L1CAM$comparison="VEH and TAD"

### REACTOME Granulocyte 1V2
##  REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX
top_gene14=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$REACTOME))
top_gene14=top_gene14[top_gene14$pathway=='REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX',]
top_gene14=unlist(top_gene14[,8])
REACTOME_DEGRAD_EXTRA=as.data.frame(Gra_1v2DE)
REACTOME_DEGRAD_EXTRA=REACTOME_DEGRAD_EXTRA[(row.names(REACTOME_DEGRAD_EXTRA) %in% top_gene14),]
REACTOME_DEGRAD_EXTRA=REACTOME_DEGRAD_EXTRA[,]
REACTOME_DEGRAD_EXTRA$pathway="REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX"
REACTOME_DEGRAD_EXTRA$celltype="Granulocyte"
REACTOME_DEGRAD_EXTRA$comparison="VEH and TAD"

### REACTOME Granulocyte 1V2
##  REACTOME_ELASTIC_FIBRE_FORMATION
top_gene15=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$REACTOME))
top_gene15=top_gene15[top_gene15$pathway=='REACTOME_ELASTIC_FIBRE_FORMATION',]
top_gene15=unlist(top_gene15[,8])
REACTOME_ELASTIC_FIBER=as.data.frame(Gra_1v2DE)
REACTOME_ELASTIC_FIBER=REACTOME_ELASTIC_FIBER[(row.names(REACTOME_ELASTIC_FIBER) %in% top_gene15),]
REACTOME_ELASTIC_FIBER=REACTOME_ELASTIC_FIBER[REACTOME_ELASTIC_FIBER$pvalue<0.05,]
REACTOME_ELASTIC_FIBER$pathway="REACTOME_ELASTIC_FIBRE_FORMATION"
REACTOME_ELASTIC_FIBER$celltype="Granulocyte"
REACTOME_ELASTIC_FIBER$comparison="VEH and TAD"

### REACTOME Granulocyte 1V2
##  REACTOME_CELL_JUNCTION_ORGANIZATION
top_gene16=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$REACTOME))
top_gene16=top_gene16[top_gene16$pathway=='REACTOME_CELL_JUNCTION_ORGANIZATION',]
top_gene16=unlist(top_gene16[,8])
REACTOME_CELL_JUNC=as.data.frame(Gra_1v2DE)
REACTOME_CELL_JUNC=REACTOME_CELL_JUNC[(row.names(REACTOME_CELL_JUNC) %in% top_gene16),]
REACTOME_CELL_JUNC=REACTOME_CELL_JUNC[REACTOME_CELL_JUNC$pvalue<0.05,]
REACTOME_CELL_JUNC$pathway="REACTOME_CELL_JUNCTION_ORGANIZATION"
REACTOME_CELL_JUNC$celltype="Granulocyte"
REACTOME_CELL_JUNC$comparison="VEH and TAD"

### REACTOME Granulocyte 1V2
##  REACTOME_CELL_CELL_JUNCTION_ORGANIZATION
top_gene17=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$REACTOME))
top_gene17=top_gene17[top_gene17$pathway=='REACTOME_CELL_CELL_JUNCTION_ORGANIZATION',]
top_gene17=unlist(top_gene17[,8])
REACTOME_CELL2_JUNC=as.data.frame(Gra_1v2DE)
REACTOME_CELL2_JUNC=REACTOME_CELL2_JUNC[(row.names(REACTOME_CELL2_JUNC) %in% top_gene17),]
REACTOME_CELL2_JUNC=REACTOME_CELL2_JUNC[REACTOME_CELL2_JUNC$pvalue<0.05,]
REACTOME_CELL2_JUNC$pathway="REACTOME_CELL_CELL_JUNCTION_ORGANIZATION"
REACTOME_CELL2_JUNC$celltype="Granulocyte"
REACTOME_CELL2_JUNC$comparison="VEH and TAD"

### GOMF Granulocyte 1v2
##  GOMF_CELL_ADHESION_MOLECULE_BINDING 
top_gene18=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$GOMF))
top_gene18=top_gene18[top_gene18$pathway=='GOMF_CELL_ADHESION_MOLECULE_BINDING',]
top_gene18=unlist(top_gene18[,8])
GOMF_CELL=as.data.frame(Gra_1v2DE)
GOMF_CELL=GOMF_CELL[(row.names(GOMF_CELL) %in% top_gene18),]
GOMF_CELL=GOMF_CELL[GOMF_CELL$pvalue<0.05,]
GOMF_CELL$pathway="GOMF_CELL_ADHESION_MOLECULE_BINDING"
GOMF_CELL$celltype="Granulocyte"
GOMF_CELL$comparison="VEH and TAD"

### NABA Granulocyte 1v2
##  NABA_SECRETED_FACTORS
top_gene20=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$NABA))
top_gene20=top_gene20[top_gene20$pathway=='NABA_SECRETED_FACTORS',]
top_gene20=unlist(top_gene20[,8])
NABA_SEC_FACTORS=as.data.frame(Gra_1v2DE)
NABA_SEC_FACTORS=NABA_SEC_FACTORS[(row.names(NABA_SEC_FACTORS) %in% top_gene20),]
NABA_SEC_FACTORS=NABA_SEC_FACTORS[NABA_SEC_FACTORS$pvalue<0.05,]
NABA_SEC_FACTORS$pathway="NABA_SECRETED_FACTORS"
NABA_SEC_FACTORS$celltype="Granulocyte"
NABA_SEC_FACTORS$comparison="VEH and TAD"

### NABA Granulocyte 1v2
##  NABA_ECM_AFFILIATED
top_gene22=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$NABA))
top_gene22=top_gene22[top_gene22$pathway=='NABA_ECM_AFFILIATED',]
top_gene22=unlist(top_gene22[,8])
NABA_ECM_AFF=as.data.frame(Gra_1v2DE)
NABA_ECM_AFF=NABA_ECM_AFF[(row.names(NABA_ECM_AFF) %in% top_gene22),]
NABA_ECM_AFF=NABA_ECM_AFF[NABA_ECM_AFF$pvalue<0.05,]
NABA_ECM_AFF$pathway="NABA_ECM_AFFILIATED"
NABA_ECM_AFF$celltype="Granulocyte"
NABA_ECM_AFF$comparison="VEH and TAD"

### NABA Granulocyte 1v2
##  NABA_MATRISOME_ASSOCIATED
top_gene23=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$NABA))
top_gene23=top_gene23[top_gene23$pathway=='NABA_MATRISOME_ASSOCIATED',]
top_gene23=unlist(top_gene23[,8])
NABA_MAT_ASS=as.data.frame(Gra_1v2DE)
NABA_MAT_ASS=NABA_MAT_ASS[(row.names(NABA_MAT_ASS) %in% top_gene23),]
NABA_MAT_ASS=NABA_MAT_ASS[NABA_MAT_ASS$pvalue<0.05,]
NABA_MAT_ASS$pathway="NABA_MATRISOME_ASSOCIATED"
NABA_MAT_ASS$celltype="Granulocyte"
NABA_MAT_ASS$comparison="VEH and TAD"

### NABA Granulocyte 1v2
##  NABA_ECM_REGULATORS
top_gene25=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$NABA))
top_gene25=top_gene25[top_gene25$pathway=='NABA_ECM_REGULATORS',]
top_gene25=unlist(top_gene25[,8])
NABA_ECM_REG=as.data.frame(Gra_1v2DE)
NABA_ECM_REG=NABA_ECM_REG[(row.names(NABA_ECM_REG) %in% top_gene25),]
NABA_ECM_REG=NABA_ECM_REG[NABA_ECM_REG$pvalue<0.05,]
NABA_ECM_REG$pathway="NABA_ECM_REGULATORS"
NABA_ECM_REG$celltype="Granulocyte"
NABA_ECM_REG$comparison="VEH and TAD"
### NABA Granulocyte 1v2
##  NABA_MATRISOME
top_gene26=as.data.frame(as.matrix(Gra_DE$ONEvTWOgs$NABA))
top_gene26=top_gene26[top_gene26$pathway=='NABA_MATRISOME',]
top_gene26=unlist(top_gene26[,8])
NABA_MAT=as.data.frame(Gra_1v2DE)
NABA_MAT=NABA_MAT[(row.names(NABA_MAT) %in% top_gene26),]
NABA_MAT=NABA_MAT[NABA_MAT$pvalue<0.05,]
NABA_MAT$pathway="NABA_MATRISOME"
NABA_MAT$celltype="Granulocyte"
NABA_MAT$comparison="VEH and TAD"

Gra1v2 <- rbind(HALL_INF_RES,
                HALL_EPI_MES_TRANS,
                HALL_IL2,
                KEGG_CYTO_REC,KEGG_LYSO,
                KEGG_LYSO,
                REACTOME_INNATE,
                REACTOME_NEUTRO,
                REACTOME_COLLAGEN_FORM,
                REACTOME_COLLAGEN_FORM,
                REACTOME_ASSEMBLY_OTHER_STRUC,
                REACTOME_EXTRA,
                REACTOME_L1CAM,
                REACTOME_DEGRAD_EXTRA,
                REACTOME_ELASTIC_FIBER,
                REACTOME_CELL_JUNC,
                REACTOME_CELL2_JUNC,
                GOMF_CELL,
                NABA_SEC_FACTORS,
                NABA_ECM_AFF,
                NABA_MAT_ASS,
                NABA_ECM_REG,
                NABA_MAT)

mean_VEH = c()
mean_TAD = c()
sub = subset(x=sobj,subset=cell_type=="Granulocyte")
sub1 = subset(x=sub,subset=group_id==c("VEH"))
sub2 = subset(x=sub,subset=group_id==c("TAD"))
for (i in 1:nrow(Gra1v2)){
    x = AverageExpression(sub1,feature=rownames(Gra1v2)[i])
    y = AverageExpression(sub2,feature=rownames(Gra1v2)[i])
    if(length(x)==0){
        mean_VEH = c(mean_VEH,"NA")
    }
    if(length(x)==0){
        mean_TAD = c(mean_TAD,"NA")
    }
    mean_VEH = c(mean_VEH,x)
    mean_TAD = c(mean_TAD,y)
}
xx = unlist(mean_VEH)
yy = unlist(mean_TAD)
g = data.frame(mean_VEH = xx,mean_TAD=yy)
Gra1v2 = cbind(Gra1v2,g)
write.csv(Gra1v2,"Gra_1v2_topGenes.csv")

##### 4v5 Granulocyte

Gra4v5 <- data.frame(baseMean=NULL,log2FoldChange=NULL,lfcSE=NULL,stat=NULL,pvalue=NULL,padj=NULL,pathway=NULL,celltype=NULL,comparison=NULL,mean_MLIST=NULL,mean_MLIST.TAD=NULL)

### HALLMARK Granulocyte 4v5
##  HALLMARK_IL2_STAT5_SIGNALING
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$HALLMARK))
top_gene=top_gene[top_gene$pathway=='HALLMARK_IL2_STAT5_SIGNALING',]
top_gene=unlist(top_gene[,8])
HALL_IL2=as.data.frame(Gra_4v5DE)
HALL_IL2=HALL_IL2[(row.names(HALL_IL2) %in% top_gene),]
HALL_IL2=HALL_IL2[HALL_IL2$pvalue<0.05,]
HALL_IL2$pathway="HALLMARK_IL2_STAT5_SIGNALING"
HALL_IL2$celltype="Granulocyte"
HALL_IL2$comparison="MLIST and MLIST.TAD"

### HALLMARK Granulocyte 4v5
##  HALLMARK_INFLAMMATORY_RESPONSE
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$HALLMARK))
top_gene=top_gene[top_gene$pathway=='HALLMARK_INFLAMMATORY_RESPONSE',]
top_gene=unlist(top_gene[,8])
HALL_INFLA_RES=as.data.frame(Gra_4v5DE)
HALL_INFLA_RES=HALL_INFLA_RES[(row.names(HALL_INFLA_RES) %in% top_gene),]
HALL_INFLA_RES=HALL_INFLA_RES[HALL_INFLA_RES$pvalue<0.05,]
HALL_INFLA_RES$pathway="HALLMARK_INFLAMMATORY_RESPONSE"
HALL_INFLA_RES$celltype="Granulocyte"
HALL_INFLA_RES$comparison="MLIST and MLIST.TAD"

### HALLMARK Granulocyte 4v5
##  HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$HALLMARK))
top_gene=top_gene[top_gene$pathway=='HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY',]
top_gene=unlist(top_gene[,8])
HALL_REA_O2=as.data.frame(Gra_4v5DE)
HALL_REA_O2=HALL_REA_O2[(row.names(HALL_REA_O2) %in% top_gene),]
HALL_REA_O2=HALL_REA_O2[HALL_REA_O2$pvalue<0.05,]
HALL_REA_O2$pathway="HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"
HALL_REA_O2$celltype="Granulocyte"
HALL_REA_O2$comparison="MLIST and MLIST.TAD"

### HALLMARK Granulocyte 4v5
##  HALLMARK_APOPTOSIS
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$HALLMARK))
top_gene=top_gene[top_gene$pathway=='HALLMARK_APOPTOSIS',]
top_gene=unlist(top_gene[,8])
HALL_APO=as.data.frame(Gra_4v5DE)
HALL_APO=HALL_APO[(row.names(HALL_APO) %in% top_gene),]
HALL_APO=HALL_APO[HALL_APO$pvalue<0.05,]
HALL_APO$pathway="HALLMARK_APOPTOSIS"
HALL_APO$celltype="Granulocyte"
HALL_APO$comparison="MLIST and MLIST.TAD"

### HALLMARK Granulocyte 4v5
##  HALLMARK_INTERFERON_ALPHA_RESPONSE
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$HALLMARK))
top_gene=top_gene[top_gene$pathway=='HALLMARK_INTERFERON_ALPHA_RESPONSE',]
top_gene=unlist(top_gene[,8])
HALL_INFI_ALPHA_RES=as.data.frame(Gra_4v5DE)
HALL_INFI_ALPHA_RES=HALL_INFI_ALPHA_RES[(row.names(HALL_INFI_ALPHA_RES) %in% top_gene),]
HALL_INFI_ALPHA_RES=HALL_INFI_ALPHA_RES[HALL_INFI_ALPHA_RES$pvalue<0.05,]
HALL_INFI_ALPHA_RES$pathway="HALLMARK_INTERFERON_ALPHA_RESPONSE"
HALL_INFI_ALPHA_RES$celltype="Granulocyte"
HALL_INFI_ALPHA_RES$comparison="MLIST and MLIST.TAD"

### HALLMARK Granulocyte 4v5
##  HALLMARK_TNFA_SIGNALING_VIA_NFKB
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$HALLMARK))
top_gene=top_gene[top_gene$pathway=='HALLMARK_TNFA_SIGNALING_VIA_NFKB',]
top_gene=unlist(top_gene[,8])
HALL_TNFA_SIG=as.data.frame(Gra_4v5DE)
HALL_TNFA_SIG=HALL_TNFA_SIG[(row.names(HALL_TNFA_SIG) %in% top_gene),]
HALL_TNFA_SIG=HALL_TNFA_SIG[HALL_TNFA_SIG$pvalue<0.05,]
HALL_TNFA_SIG$pathway="HALLMARK_TNFA_SIGNALING_VIA_NFKB"
HALL_TNFA_SIG$celltype="Granulocyte"
HALL_TNFA_SIG$comparison="MLIST and MLIST.TAD"

### HALLMARK Granulocyte 4v5
##  HALLMARK_INTERFERON_GAMMA_RESPONSE
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$HALLMARK))
top_gene=top_gene[top_gene$pathway=='HALLMARK_INTERFERON_GAMMA_RESPONSE',]
top_gene=unlist(top_gene[,8])
HALL_INFI_GAMMA_RES=as.data.frame(Gra_4v5DE)
HALL_INFI_GAMMA_RES=HALL_INFI_GAMMA_RES[(row.names(HALL_INFI_GAMMA_RES) %in% top_gene),]
HALL_INFI_GAMMA_RES=HALL_INFI_GAMMA_RES[HALL_INFI_GAMMA_RES$pvalue<0.05,]
HALL_INFI_GAMMA_RES$pathway="HALLMARK_INTERFERON_GAMMA_RESPONSE"
HALL_INFI_GAMMA_RES$celltype="Granulocyte"
HALL_INFI_GAMMA_RES$comparison="MLIST and MLIST.TAD"

### HALLMARK Granulocyte 4v5
##  HALLMARK_HYPOXIA
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$HALLMARK))
top_gene=top_gene[top_gene$pathway=='HALLMARK_HYPOXIA',]
top_gene=unlist(top_gene[,8])
HALL_HYPOXIA=as.data.frame(Gra_4v5DE)
HALL_HYPOXIA=HALL_HYPOXIA[(row.names(HALL_HYPOXIA) %in% top_gene),]
HALL_HYPOXIA=HALL_HYPOXIA[HALL_HYPOXIA$pvalue<0.05,]
HALL_HYPOXIA$pathway="HALLMARK_HYPOXIA"
HALL_HYPOXIA$celltype="Granulocyte"
HALL_HYPOXIA$comparison="MLIST and MLIST.TAD"

### KEGG Granulocyte 4V5
##  KEGG_PROTEASOME
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$KEGG))
top_gene=top_gene[top_gene$pathway=='KEGG_PROTEASOME',]
top_gene=unlist(top_gene[,8])
KEGG_PRO=as.data.frame(Gra_4v5DE)
KEGG_PRO=KEGG_PRO[(row.names(KEGG_PRO) %in% top_gene),]
KEGG_PRO=KEGG_PRO[KEGG_PRO$pvalue<0.05,]
KEGG_PRO$pathway="KEGG_PROTEASOME"
KEGG_PRO$celltype="Granulocyte"
KEGG_PRO$comparison="MLIST and MLIST.TAD"

### KEGG Granulocyte 4V5
##  KEGG_LYSOSOME
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$KEGG))
top_gene=top_gene[top_gene$pathway=='KEGG_LYSOSOME',]
top_gene=unlist(top_gene[,8])
KEGG_LYSO=as.data.frame(Gra_4v5DE)
KEGG_LYSO=KEGG_LYSO[(row.names(KEGG_LYSO) %in% top_gene),]
KEGG_LYSO=KEGG_LYSO[KEGG_LYSO$pvalue<0.05,]
KEGG_LYSO$pathway="KEGG_LYSOSOME"
KEGG_LYSO$celltype="Granulocyte"
KEGG_LYSO$comparison="MLIST and MLIST.TAD"

### KEGG Granulocyte 4V5
##  KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$KEGG))
top_gene=top_gene[top_gene$pathway=='KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY',]
top_gene=unlist(top_gene[,8])
KEGG_NOD_LIKE=as.data.frame(Gra_4v5DE)
KEGG_NOD_LIKE=KEGG_NOD_LIKE[(row.names(KEGG_NOD_LIKE) %in% top_gene),]
KEGG_NOD_LIKE=KEGG_NOD_LIKE[KEGG_NOD_LIKE$pvalue<0.05,]
KEGG_NOD_LIKE$pathway="KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY"
KEGG_NOD_LIKE$celltype="Granulocyte"
KEGG_NOD_LIKE$comparison="MLIST and MLIST.TAD"

### KEGG Granulocyte 4V5
##  KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$KEGG))
top_gene=top_gene[top_gene$pathway=='KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY',]
top_gene=unlist(top_gene[,8])
KEGG_TOLL_LIKE=as.data.frame(Gra_4v5DE)
KEGG_TOLL_LIKE=KEGG_TOLL_LIKE[(row.names(KEGG_TOLL_LIKE) %in% top_gene),]
KEGG_TOLL_LIKE=KEGG_TOLL_LIKE[KEGG_TOLL_LIKE$pvalue<0.05,]
KEGG_TOLL_LIKE$pathway="KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY"
KEGG_TOLL_LIKE$celltype="Granulocyte"
KEGG_TOLL_LIKE$comparison="MLIST and MLIST.TAD"

### KEGG Granulocyte 4V5
##  KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$KEGG))
top_gene=top_gene[top_gene$pathway=='KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY',]
top_gene=unlist(top_gene[,8])
KEGG_RIG_I=as.data.frame(Gra_4v5DE)
KEGG_RIG_I=KEGG_RIG_I[(row.names(KEGG_RIG_I) %in% top_gene),]
KEGG_RIG_I=KEGG_RIG_I[KEGG_RIG_I$pvalue<0.05,]
KEGG_RIG_I$pathway="KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY"
KEGG_RIG_I$celltype="Granulocyte"
KEGG_RIG_I$comparison="MLIST and MLIST.TAD"

### REACTOME Granulocyte 4V5
##  REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$REACTOME))
top_gene=top_gene[top_gene$pathway=='REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM',]
top_gene=unlist(top_gene[,8])
REACTOME_CYTOK=as.data.frame(Gra_4v5DE)
REACTOME_CYTOK=REACTOME_CYTOK[(row.names(REACTOME_CYTOK) %in% top_gene),]
REACTOME_CYTOK=REACTOME_CYTOK[REACTOME_CYTOK$pvalue<0.05,]
REACTOME_CYTOK$pathway="REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM"
REACTOME_CYTOK$celltype="Granulocyte"
REACTOME_CYTOK$comparison="MLIST and MLIST.TAD"

### REACTOME Granulocyte 4V5
##  REACTOME_INTERFERON_SIGNALING
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$REACTOME))
top_gene=top_gene[top_gene$pathway=='REACTOME_INTERFERON_SIGNALING',]
top_gene=unlist(top_gene[,8])
REACTOME_INTER=as.data.frame(Gra_4v5DE)
REACTOME_INTER=REACTOME_INTER[(row.names(REACTOME_INTER) %in% top_gene),]
REACTOME_INTER=REACTOME_INTER[REACTOME_INTER$pvalue<0.05,]
REACTOME_INTER$pathway="REACTOME_INTERFERON_SIGNALING"
REACTOME_INTER$celltype="Granulocyte"
REACTOME_INTER$comparison="MLIST and MLIST.TAD"

### REACTOME Granulocyte 4V5
##  REACTOME_TOLL_LIKE_RECEPTOR_CASCADES
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$REACTOME))
top_gene=top_gene[top_gene$pathway=='REACTOME_TOLL_LIKE_RECEPTOR_CASCADES',]
top_gene=unlist(top_gene[,8])
REACTOME_TOLLLIKE=as.data.frame(Gra_4v5DE)
REACTOME_TOLLLIKE=REACTOME_TOLLLIKE[(row.names(REACTOME_TOLLLIKE) %in% top_gene),]
REACTOME_TOLLLIKE=REACTOME_TOLLLIKE[REACTOME_TOLLLIKE$pvalue<0.05,]
REACTOME_TOLLLIKE$pathway="REACTOME_TOLL_LIKE_RECEPTOR_CASCADES"
REACTOME_TOLLLIKE$celltype="Granulocyte"
REACTOME_TOLLLIKE$comparison="MLIST and MLIST.TAD"

### REACTOME Granulocyte 4V5
##  REACTOME_NEUTROPHIL_DEGRANULATION
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$REACTOME))
top_gene=top_gene[top_gene$pathway=='REACTOME_NEUTROPHIL_DEGRANULATION',]
top_gene=unlist(top_gene[,8])
REACTOME_NEWTRO=as.data.frame(Gra_4v5DE)
REACTOME_NEWTRO=REACTOME_NEWTRO[(row.names(REACTOME_NEWTRO) %in% top_gene),]
REACTOME_NEWTRO=REACTOME_NEWTRO[REACTOME_NEWTRO$pvalue<0.05,]
REACTOME_NEWTRO$pathway="REACTOME_NEUTROPHIL_DEGRANULATION"
REACTOME_NEWTRO$celltype="Granulocyte"
REACTOME_NEWTRO$comparison="MLIST and MLIST.TAD"

### REACTOME Granulocyte 4V5
##  REACTOME_INNATE_IMMUNE_SYSTEM
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$REACTOME))
top_gene=top_gene[top_gene$pathway=='REACTOME_INNATE_IMMUNE_SYSTEM',]
top_gene=unlist(top_gene[,8])
REACTOME_INNATE=as.data.frame(Gra_4v5DE)
REACTOME_INNATE=REACTOME_INNATE[(row.names(REACTOME_INNATE) %in% top_gene),]
REACTOME_INNATE=REACTOME_INNATE[REACTOME_INNATE$pvalue<0.05,]
REACTOME_INNATE$pathway="REACTOME_INNATE_IMMUNE_SYSTEM"
REACTOME_INNATE$celltype="Granulocyte"
REACTOME_INNATE$comparison="MLIST and MLIST.TAD"

### REACTOME Granulocyte 4V5
##  REACTOME_MYD88_INDEPENDENT_TLR4_CASCADE
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$REACTOME))
top_gene=top_gene[top_gene$pathway=='REACTOME_MYD88_INDEPENDENT_TLR4_CASCADE',]
top_gene=unlist(top_gene[,8])
REACTOME_MYD88=as.data.frame(Gra_4v5DE)
REACTOME_MYD88=REACTOME_MYD88[(row.names(REACTOME_MYD88) %in% top_gene),]
REACTOME_MYD88=REACTOME_MYD88[REACTOME_MYD88$pvalue<0.05,]
REACTOME_MYD88$pathway="REACTOME_MYD88_INDEPENDENT_TLR4_CASCADE"
REACTOME_MYD88$celltype="Granulocyte"
REACTOME_MYD88$comparison="MLIST and MLIST.TAD"

### REACTOME Granulocyte 4V5
##  REACTOME_ROS_AND_RNS_PRODUCTION_IN_PHAGOCYTES
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$REACTOME))
top_gene=top_gene[top_gene$pathway=='REACTOME_ROS_AND_RNS_PRODUCTION_IN_PHAGOCYTES',]
top_gene=unlist(top_gene[,8])
REACTOME_ROS_RNS=as.data.frame(Gra_4v5DE)
REACTOME_ROS_RNS=REACTOME_ROS_RNS[(row.names(REACTOME_ROS_RNS) %in% top_gene),]
REACTOME_ROS_RNS=REACTOME_ROS_RNS[REACTOME_ROS_RNS$pvalue<0.05,]
REACTOME_ROS_RNS$pathway="REACTOME_ROS_AND_RNS_PRODUCTION_IN_PHAGOCYTES"
REACTOME_ROS_RNS$celltype="Granulocyte"
REACTOME_ROS_RNS$comparison="MLIST and MLIST.TAD"

### REACTOME Granulocyte 4V5
##  REACTOME_THE_NLRP3_INFLAMMASOME
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$REACTOME))
top_gene=top_gene[top_gene$pathway=='REACTOME_THE_NLRP3_INFLAMMASOME',]
top_gene=unlist(top_gene[,8])
REACTOME_NLRP3=as.data.frame(Gra_4v5DE)
REACTOME_NLRP3=REACTOME_NLRP3[(row.names(REACTOME_NLRP3) %in% top_gene),]
REACTOME_NLRP3=REACTOME_NLRP3[REACTOME_NLRP3$pvalue<0.05,]
REACTOME_NLRP3$pathway="REACTOME_THE_NLRP3_INFLAMMASOME"
REACTOME_NLRP3$celltype="Granulocyte"
REACTOME_NLRP3$comparison="MLIST and MLIST.TAD"

### REACTOME Granulocyte 4V5
##  REACTOME_TRAF3_DEPENDENT_IRF_ACTIVATION_PATHWAY
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$REACTOME))
top_gene=top_gene[top_gene$pathway=='REACTOME_TRAF3_DEPENDENT_IRF_ACTIVATION_PATHWAY',]
top_gene=unlist(top_gene[,8])
REACTOME_TRAF3=as.data.frame(Gra_4v5DE)
REACTOME_TRAF3=REACTOME_TRAF3[(row.names(REACTOME_TRAF3) %in% top_gene),]
REACTOME_TRAF3=REACTOME_TRAF3[REACTOME_TRAF3$pvalue<0.05,]
REACTOME_TRAF3$pathway="REACTOME_TRAF3_DEPENDENT_IRF_ACTIVATION_PATHWAY"
REACTOME_TRAF3$celltype="Granulocyte"
REACTOME_TRAF3$comparison="MLIST and MLIST.TAD"

### REACTOME Granulocyte 4V5
##  REACTOME_OAS_ANTIVIRAL_RESPONSE
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$REACTOME))
top_gene=top_gene[top_gene$pathway=='REACTOME_OAS_ANTIVIRAL_RESPONSE',]
top_gene=unlist(top_gene[,8])
REACTOME_OAS=as.data.frame(Gra_4v5DE)
REACTOME_OAS=REACTOME_OAS[(row.names(REACTOME_OAS) %in% top_gene),]
REACTOME_OAS=REACTOME_OAS[REACTOME_OAS$pvalue<0.05,]
REACTOME_OAS$pathway="REACTOME_OAS_ANTIVIRAL_RESPONSE"
REACTOME_OAS$celltype="Granulocyte"
REACTOME_OAS$comparison="MLIST and MLIST.TAD"

### REACTOME Granulocyte 4V5
##  REACTOME_DDX58_IFIH1_MEDIATED_INDUCTION_OF_INTERFERON_ALPHA_BETA
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$REACTOME))
top_gene=top_gene[top_gene$pathway=='REACTOME_DDX58_IFIH1_MEDIATED_INDUCTION_OF_INTERFERON_ALPHA_BETA',]
top_gene=unlist(top_gene[,8])
REACTOME_DDX58=as.data.frame(Gra_4v5DE)
REACTOME_DDX58=REACTOME_DDX58[(row.names(REACTOME_DDX58) %in% top_gene),]
REACTOME_DDX58=REACTOME_DDX58[REACTOME_DDX58$pvalue<0.05,]
REACTOME_DDX58$pathway="REACTOME_DDX58_IFIH1_MEDIATED_INDUCTION_OF_INTERFERON_ALPHA_BETA"
REACTOME_DDX58$celltype="Granulocyte"
REACTOME_DDX58$comparison="MLIST and MLIST.TAD"

### REACTOME Granulocyte 4V5
##  REACTOME_INFLAMMASOMES
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$REACTOME))
top_gene=top_gene[top_gene$pathway=='REACTOME_INFLAMMASOMES',]
top_gene=unlist(top_gene[,8])
REACTOME_INFLA=as.data.frame(Gra_4v5DE)
REACTOME_INFLA=REACTOME_INFLA[(row.names(REACTOME_INFLA) %in% top_gene),]
REACTOME_INFLA=REACTOME_INFLA[REACTOME_INFLA$pvalue<0.05,]
REACTOME_INFLA$pathway="REACTOME_INFLAMMASOMES"
REACTOME_INFLA$celltype="Granulocyte"
REACTOME_INFLA$comparison="MLIST and MLIST.TAD"

### REACTOME Granulocyte 4V5
##  REACTOME_INTERFERON_ALPHA_BETA_SIGNALING
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$REACTOME))
top_gene=top_gene[top_gene$pathway=='REACTOME_INTERFERON_ALPHA_BETA_SIGNALING',]
top_gene=unlist(top_gene[,8])
REACTOME_INTER_ALPHA=as.data.frame(Gra_4v5DE)
REACTOME_INTER_ALPHA=REACTOME_INTER_ALPHA[(row.names(REACTOME_INTER_ALPHA) %in% top_gene),]
REACTOME_INTER_ALPHA=REACTOME_INTER_ALPHA[REACTOME_INTER_ALPHA$pvalue<0.05,]
REACTOME_INTER_ALPHA$pathway="REACTOME_INTERFERON_ALPHA_BETA_SIGNALING"
REACTOME_INTER_ALPHA$celltype="Granulocyte"
REACTOME_INTER_ALPHA$comparison="MLIST and MLIST.TAD"

### NABA Granulocyte 4v5
##  NABA_ECM_REGULATORS
top_gene=as.data.frame(as.matrix(Gra_DE$FOURvFIVEgs$NABA))
top_gene=top_gene[top_gene$pathway=='NABA_ECM_REGULATORS',]
top_gene=unlist(top_gene[,8])
NABA_ECM=as.data.frame(Gra_4v5DE)
NABA_ECM=NABA_ECM[(row.names(NABA_ECM) %in% top_gene),]
NABA_ECM=NABA_ECM[NABA_ECM$pvalue<0.05,]
NABA_ECM$pathway="NABA_ECM_REGULATORS"
NABA_ECM$celltype="Granulocyte"
NABA_ECM$comparison="MLIST and MLIST.TAD"

Gra4v5 <- rbind(HALL_IL2,
                HALL_INFLA_RES,
                HALL_REA_O2,
                HALL_APO,
                HALL_INFI_ALPHA_RES,
                HALL_TNFA_SIG,
                HALL_INFI_GAMMA_RES,
                HALL_HYPOXIA,
                KEGG_PRO,
                KEGG_LYSO,
                KEGG_NOD_LIKE,
                KEGG_TOLL_LIKE,
                KEGG_RIG_I,
                REACTOME_CYTOK,
                REACTOME_INTER,
                REACTOME_TOLLLIKE,
                REACTOME_NEWTRO,
                REACTOME_INNATE,
                REACTOME_MYD88,
                REACTOME_ROS_RNS,
                REACTOME_NLRP3,
                REACTOME_TRAF3,
                REACTOME_OAS,
                REACTOME_DDX58,
                REACTOME_INFLA,
                REACTOME_INTER_ALPHA,
                NABA_ECM
               )

mean_MLIST = c()
mean_MLIST.TAD = c()
sub = subset(x=sobj,subset=cell_type=="Granulocyte")
sub1 = subset(x=sub,subset=group_id==c("MLIST"))
sub2 = subset(x=sub,subset=group_id==c("MLIST.TAD"))
for (i in 1:nrow(Gra4v5)){
    x = AverageExpression(sub1,feature=rownames(Gra4v5)[i])
    y = AverageExpression(sub2,feature=rownames(Gra4v5)[i])
    if(length(x)==0){
        mean_MLIST = c(mean_MLIST,"NA")
    }
    if(length(x)==0){
        mean_MLIST.TAD = c(mean_MLIST.TAD,"NA")
    }
    mean_MLIST = c(mean_MLIST,x)
    mean_MLIST.TAD = c(mean_MLIST.TAD,y)
}
xx = unlist(mean_MLIST)
yy = unlist(mean_MLIST.TAD)
g = data.frame(mean_MLIST = xx,mean_MLIST.TAD=yy)
Gra4v5 = cbind(Gra4v5,g)
write.csv(Gra4v5,"Gra_4v5_topGenes.csv")

##### 1v2 Macrophage

Mac1v2 <- data.frame(baseMean=NULL,log2FoldChange=NULL,lfcSE=NULL,stat=NULL,pvalue=NULL,padj=NULL,pathway=NULL,celltype=NULL,comparison=NULL,mean_VEH=NULL,mean_TAD=NULL)

### HALLMARK Macrophage 1v2
##  HALLMARK_IL2_STAT5_SIGNALING
top_gene=as.data.frame(as.matrix(Mac_DE$ONEvTWOgs$HALLMARK))
top_gene=top_gene[top_gene$pathway=='HALLMARK_IL2_STAT5_SIGNALING',]
top_gene=unlist(top_gene[,8])
HALL_IL2_STAT5_SIG=as.data.frame(Mac_1v2DE)
HALL_IL2_STAT5_SIG=HALL_IL2_STAT5_SIG[(row.names(HALL_IL2_STAT5_SIG) %in% top_gene),]
HALL_IL2_STAT5_SIG=HALL_IL2_STAT5_SIG[HALL_IL2_STAT5_SIG$pvalue<0.05,]
HALL_IL2_STAT5_SIG$pathway="HALLMARK_IL2_STAT5_SIGNALING"
HALL_IL2_STAT5_SIG$celltype="Macrophage"
HALL_IL2_STAT5_SIG$comparison="VEH and TAD"

### NABA Macrophage 1v2
##  NABA_MATRISOME
top_gene=as.data.frame(as.matrix(Mac_DE$ONEvTWOgs$NABA))
top_gene=top_gene[top_gene$pathway=='NABA_MATRISOME',]
top_gene=unlist(top_gene[,8])
NABA_ECM=as.data.frame(Mac_1v2DE)
NABA_ECM=NABA_ECM[(row.names(NABA_ECM) %in% top_gene),]
NABA_ECM=NABA_ECM[NABA_ECM$pvalue<0.05,]


### NABA Macrophage 1v2
##  NABA_MATRISOME_ASSOCIATED
top_gene=as.data.frame(as.matrix(Mac_DE$ONEvTWOgs$NABA))
top_gene=top_gene[top_gene$pathway=='NABA_MATRISOME_ASSOCIATED',]
top_gene=unlist(top_gene[,8])
NABA_MATRI=as.data.frame(Mac_1v2DE)
NABA_MATRI=NABA_MATRI[(row.names(NABA_MATRI) %in% top_gene),]
NABA_MATRI=NABA_MATRI[NABA_MATRI$pvalue<0.05,]


### NABA Macrophage 1v2
##  NABA_SECRETED_FACTORS
top_gene=as.data.frame(as.matrix(Mac_DE$ONEvTWOgs$NABA))
top_gene=top_gene[top_gene$pathway=='NABA_SECRETED_FACTORS',]
top_gene=unlist(top_gene[,8])
NABA_SEC_FACTOR=as.data.frame(Mac_1v2DE)
NABA_SEC_FACTOR=NABA_SEC_FACTOR[(row.names(NABA_SEC_FACTOR) %in% top_gene),]
NABA_SEC_FACTOR=NABA_SEC_FACTOR[NABA_SEC_FACTOR$pvalue<0.05,]

Mac1v2 <- rbind(HALL_IL2_STAT5_SIG,
                NABA_MATRI,
                NABA_SEC_FACTOR,
                NABA_ECM
               ) 

mean_VEH = c()
mean_TAD = c()
sub = subset(x=sobj,subset=cell_type=="Macrophage")
sub1 = subset(x=sub,subset=group_id==c("VEH"))
sub2 = subset(x=sub,subset=group_id==c("TAD"))
for (i in 1:nrow(Mac1v2)){
    x = AverageExpression(sub1,feature=rownames(Mac1v2)[i])
    y = AverageExpression(sub2,feature=rownames(Mac1v2)[i])
    if(length(x)==0){
        mean_VEH = c(mean_VEH,"NA")
    }
    if(length(x)==0){
        mean_TAD = c(mean_TAD,"NA")
    }
    mean_VEH = c(mean_VEH,x)
    mean_TAD = c(mean_TAD,y)
}
xx = unlist(mean_VEH)
yy = unlist(mean_TAD)
g = data.frame(mean_VEH = xx,mean_TAD=yy)
Mac1v2 = cbind(Mac1v2,g)
write.csv(Mac1v2,"Mac_1v2_topGenes.csv")

##### 4v5 Macrophage

Mac4v5 <- data.frame(baseMean=NULL,log2FoldChange=NULL,lfcSE=NULL,stat=NULL,pvalue=NULL,padj=NULL,pathway=NULL,celltype=NULL,comparison=NULL,mean_MLIST=NULL,mean_MLIST.TAD=NULL)

### HALLMARK Macrophage 4v5
##  HALLMARK_INTERFERON_ALPHA_RESPONSE
top_gene=as.data.frame(as.matrix(Mac_DE$FOURvFIVEgs$HALLMARK))
top_gene=top_gene[top_gene$pathway=='HALLMARK_INTERFERON_ALPHA_RESPONSE',]
top_gene=unlist(top_gene[,8])
HALL_INTER_ALPHA_RES=as.data.frame(Mac_4v5DE)
HALL_INTER_ALPHA_RES=HALL_INTER_ALPHA_RES[(row.names(HALL_INTER_ALPHA_RES) %in% top_gene),]
HALL_INTER_ALPHA_RES=HALL_INTER_ALPHA_RES[HALL_INTER_ALPHA_RES$pvalue<0.05,]
if (nrow(HALL_INTER_ALPHA_RES)!=0){
    HALL_INTER_ALPHA_RES$pathway="HALLMARK_INTERFERON_ALPHA_RESPONSE"
    HALL_INTER_ALPHA_RES$celltype="Macrophage"
    HALL_INTER_ALPHA_RES$comparison="MLIST and MLIST.TAD"
}
### HALLMARK Macrophage 4v5
##  HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY
top_gene=as.data.frame(as.matrix(Mac_DE$FOURvFIVEgs$HALLMARK))
top_gene=top_gene[top_gene$pathway=='HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY',]
top_gene=unlist(top_gene[,8])
HALL_REAC_O2=as.data.frame(Mac_4v5DE)
HALL_REAC_O2=HALL_REAC_O2[(row.names(HALL_REAC_O2) %in% top_gene),]
HALL_REAC_O2=HALL_REAC_O2[HALL_REAC_O2$pvalue<0.05,]
if (nrow(HALL_REAC_O2)!=0){
    HALL_REAC_O2$pathway="HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"
    HALL_REAC_O2$celltype="Macrophage"
    HALL_REAC_O2$comparison="MLIST and MLIST.TAD"
}
### HALLMARK Macrophage 4v5
##  HALLMARK_HYPOXIA
top_gene=as.data.frame(as.matrix(Mac_DE$FOURvFIVEgs$HALLMARK))
top_gene=top_gene[top_gene$pathway=='HALLMARK_HYPOXIA',]
top_gene=unlist(top_gene[,8])
HALL_HYPOXIA=as.data.frame(Mac_4v5DE)
HALL_HYPOXIA=HALL_HYPOXIA[(row.names(HALL_HYPOXIA) %in% top_gene),]
HALL_HYPOXIA=HALL_HYPOXIA[HALL_HYPOXIA$pvalue<0.05,]
if (nrow(HALL_HYPOXIA)!=0){
    HALL_HYPOXIA$pathway="HALLMARK_HYPOXIA"
    HALL_HYPOXIA$celltype="Macrophage"
    HALL_HYPOXIA$comparison="MLIST and MLIST.TAD"
}
### HALLMARK Macrophage 4v5
##  HALLMARK_INTERFERON_GAMMA_RESPONSE
top_gene=as.data.frame(as.matrix(Mac_DE$FOURvFIVEgs$HALLMARK))
top_gene=top_gene[top_gene$pathway=='HALLMARK_INTERFERON_GAMMA_RESPONSE',]
top_gene=unlist(top_gene[,8])
HALL_INTER_GAMMA_RES=as.data.frame(Mac_4v5DE)
HALL_INTER_GAMMA_RES=HALL_INTER_GAMMA_RES[(row.names(HALL_INTER_GAMMA_RES) %in% top_gene),]
HALL_INTER_GAMMA_RES=HALL_INTER_GAMMA_RES[HALL_INTER_GAMMA_RES$pvalue<0.05,]
if (nrow(HALL_INTER_GAMMA_RES)!=0){
    HALL_INTER_GAMMA_RES$pathway="HALLMARK_INTERFERON_GAMMA_RESPONSE"
    HALL_INTER_GAMMA_RES$celltype="Macrophage"
    HALL_INTER_GAMMA_RES$comparison="MLIST and MLIST.TAD"
}
### HALLMARK Macrophage 4V5
##  HALLMARK_IL2_STAT5_SIGNALING
top_gene=as.data.frame(as.matrix(Mac_DE$FOURvFIVEgs$HALLMARK))
top_gene=top_gene[top_gene$pathway=='HALLMARK_IL2_STAT5_SIGNALING',]
top_gene=unlist(top_gene[,8])
HALL_IL2_STAT5_SIG=as.data.frame(Mac_4v5DE)
HALL_IL2_STAT5_SIG=HALL_IL2_STAT5_SIG[(row.names(HALL_IL2_STAT5_SIG) %in% top_gene),]
HALL_IL2_STAT5_SIG=HALL_IL2_STAT5_SIG[HALL_IL2_STAT5_SIG$pvalue<0.05,]
if (nrow(HALL_IL2_STAT5_SIG)!=0){
    HALL_IL2_STAT5_SIG$pathway="HALLMARK_IL2_STAT5_SIGNALING"
    HALL_IL2_STAT5_SIG$celltype="Macrophage"
    HALL_IL2_STAT5_SIG$comparison="MLIST and MLIST.TAD"
}
### KEGG Macrophage 4V5
##  KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION
top_gene=as.data.frame(as.matrix(Mac_DE$FOURvFIVEgs$KEGG))
top_gene=top_gene[top_gene$pathway=='KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION',]
top_gene=unlist(top_gene[,8])
KEGG_ANTI_PRO=as.data.frame(Mac_4v5DE)
KEGG_ANTI_PRO=KEGG_ANTI_PRO[(row.names(KEGG_ANTI_PRO) %in% top_gene),]
KEGG_ANTI_PRO=KEGG_ANTI_PRO[KEGG_ANTI_PRO$pvalue<0.05,]
if (nrow(KEGG_ANTI_PRO)!=0){
    KEGG_ANTI_PRO$pathway="KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION"
    KEGG_ANTI_PRO$celltype="Macrophage"
    KEGG_ANTI_PRO$comparison="MLIST and MLIST.TAD"
}
Mac4v5 <- rbind(HALL_INTER_ALPHA_RES,
                HALL_REAC_O2,
                HALL_HYPOXIA,
                HALL_INTER_GAMMA_RES,
                HALL_IL2_STAT5_SIG,
                KEGG_ANTI_PRO
               )

mean_MLIST = c()
mean_MLIST.TAD = c()
sub = subset(x=sobj,subset=cell_type=="Granulocyte")
sub1 = subset(x=sub,subset=group_id==c("MLIST"))
sub2 = subset(x=sub,subset=group_id==c("MLIST.TAD"))
for (i in 1:nrow(Gra4v5)){
    x = AverageExpression(sub1,feature=rownames(Gra4v5)[i])
    y = AverageExpression(sub2,feature=rownames(Gra4v5)[i])
    if(length(x)==0){
        mean_MLIST = c(mean_MLIST,"NA")
    }
    if(length(x)==0){
        mean_MLIST.TAD = c(mean_MLIST.TAD,"NA")
    }
    mean_MLIST = c(mean_MLIST,x)
    mean_MLIST.TAD = c(mean_MLIST.TAD,y)
}
xx = unlist(mean_MLIST)
yy = unlist(mean_MLIST.TAD)
g = data.frame(mean_MLIST = xx,mean_MLIST.TAD=yy)
Gra4v5 = cbind(Gra4v5,g)
write.csv(Gra4v5,"Gra_4v5_topGenes.csv")

##### 1v2 Tcell

Tcell1v2 <- data.frame(baseMean=NULL,log2FoldChange=NULL,lfcSE=NULL,stat=NULL,pvalue=NULL,padj=NULL,pathway=NULL,celltype=NULL,comparison=NULL)

### HALLMARK Tcell 1v2
##  HALLMARK_INFLAMMATORY_RESPONSE
top_gene=as.data.frame(as.matrix(Tcell_DE$ONEvTWOgs$HALLMARK))
top_gene=top_gene[top_gene$pathway=='HALLMARK_INFLAMMATORY_RESPONSE',]
top_gene=unlist(top_gene[,8])
HALL_INFLA_RES=as.data.frame(Tcell_1v2DE)
HALL_INFLA_RES=HALL_INFLA_RES[(row.names(HALL_INFLA_RES) %in% top_gene),]
HALL_INFLA_RES=HALL_INFLA_RES[HALL_INFLA_RES$pvalue<0.05,]
if (nrow(HALL_INFLA_RES)!=0){
    HALL_INFLA_RES$pathway="HALLMARK_INFLAMMATORY_RESPONSE"
    HALL_INFLA_RES$celltype="Tcell"
    HALL_INFLA_RES$comparison="VEH and TAD"
}
### GOMF Tcell 1v2
##  GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT
top_gene=as.data.frame(as.matrix(Tcell_DE$ONEvTWOgs$GOMF))
top_gene=top_gene[top_gene$pathway=='GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT',]
top_gene=unlist(top_gene[,8])
GOMF_EXTRA_tcell=as.data.frame(Tcell_1v2DE)
GOMF_EXTRA_tcell=GOMF_EXTRA_tcell[(row.names(GOMF_EXTRA_tcell) %in% top_gene),]
GOMF_EXTRA_tcell=GOMF_EXTRA_tcell[GOMF_EXTRA_tcell$pvalue<0.05,]
if (nrow(GOMF_EXTRA_tcell)!=0){
    GOMF_EXTRA_tcell$pathway="GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT"
    GOMF_EXTRA_tcell$celltype="Tcell"
    GOMF_EXTRA_tcell$comparison="VEH and TAD"
}
### NABA Tcell 1v2
##  NABA_ECM_REGULATORS
top_gene=as.data.frame(as.matrix(Tcell_DE$ONEvTWOgs$NABA))
top_gene=top_gene[top_gene$pathway=='NABA_ECM_REGULATORS',]
top_gene=unlist(top_gene[,8])
NABA_ECM_REG=as.data.frame(Tcell_1v2DE)
NABA_ECM_REG=NABA_ECM_REG[(row.names(NABA_ECM_REG) %in% top_gene),]
NABA_ECM_REG=NABA_ECM_REG[NABA_ECM_REG$pvalue<0.05,]
if (nrow(NABA_ECM_REG)!=0){
    NABA_ECM_REG$pathway="NABA_ECM_REGULATORS"
    NABA_ECM_REG$celltype="Tcell"
    NABA_ECM_REG$comparison="VEH and TAD"
}
### NABA Tcell 1v2
##  NABA_SECRETED_FACTORS
top_gene=as.data.frame(as.matrix(Tcell_DE$ONEvTWOgs$NABA))
top_gene=top_gene[top_gene$pathway=='NABA_SECRETED_FACTORS',]
top_gene=unlist(top_gene[,8])
NABA_SECRET_FAC=as.data.frame(Tcell_1v2DE)
NABA_SECRET_FAC=NABA_SECRET_FAC[(row.names(NABA_SECRET_FAC) %in% top_gene),]
NABA_SECRET_FAC=NABA_SECRET_FAC[NABA_SECRET_FAC$pvalue<0.05,]
if (nrow(NABA_SECRET_FAC)!=0){
    NABA_SECRET_FAC$pathway="NABA_SECRETED_FACTORS"
    NABA_SECRET_FAC$celltype="Tcell"
    NABA_SECRET_FAC$comparison="VEH and TAD"
}
### NABA Tcell 1v2
##  NABA_ECM_AFFILIATED
top_gene=as.data.frame(as.matrix(Tcell_DE$ONEvTWOgs$NABA))
top_gene=top_gene[top_gene$pathway=='NABA_ECM_AFFILIATED',]
top_gene=unlist(top_gene[,8])
NABA_ECM_AFF=as.data.frame(Tcell_1v2DE)
NABA_ECM_AFF=NABA_ECM_AFF[(row.names(NABA_ECM_AFF) %in% top_gene),]
NABA_ECM_AFF=NABA_ECM_AFF[NABA_ECM_AFF$pvalue<0.05,]
if (nrow(NABA_ECM_AFF)!=0){
    NABA_ECM_AFF$pathway="NABA_ECM_AFFILIATED"
    NABA_ECM_AFF$celltype="Tcell"
    NABA_ECM_AFF$comparison="VEH and TAD"
}
### NABA Tcell 1v2
##  NABA_COLLAGENS
top_gene=as.data.frame(as.matrix(Tcell_DE$ONEvTWOgs$NABA))
top_gene=top_gene[top_gene$pathway=='NABA_COLLAGENS',]
top_gene=unlist(top_gene[,8])
NABA_COLL=as.data.frame(Tcell_1v2DE)
NABA_COLL=NABA_COLL[(row.names(NABA_COLL) %in% top_gene),]
NABA_COLL=NABA_COLL[NABA_COLL$pvalue<0.05,]
if (nrow(NABA_COLL)!=0){
    NABA_COLL$pathway="NABA_COLLAGENS"
    NABA_COLL$celltype="Tcell"
    NABA_COLL$comparison="VEH and TAD"
}
### NABA Tcell 1v2
##  NABA_CORE_MATRISOME
top_gene=as.data.frame(as.matrix(Tcell_DE$ONEvTWOgs$NABA))
top_gene=top_gene[top_gene$pathway=='NABA_CORE_MATRISOME',]
top_gene=unlist(top_gene[,8])
NABA_CORE_MAT=as.data.frame(Tcell_1v2DE)
NABA_CORE_MAT=NABA_CORE_MAT[(row.names(NABA_CORE_MAT) %in% top_gene),]
NABA_CORE_MAT=NABA_CORE_MAT[NABA_CORE_MAT$pvalue<0.05,]
if (nrow(NABA_CORE_MAT)!=0){
    NABA_CORE_MAT$pathway="NABA_CORE_MATRISOME"
    NABA_CORE_MAT$celltype="Tcell"
    NABA_CORE_MAT$comparison="VEH and TAD"
}
### NABA Tcell 1v2
##  NABA_MATRISOME_ASSOCIATED
top_gene=as.data.frame(as.matrix(Tcell_DE$ONEvTWOgs$NABA))
top_gene=top_gene[top_gene$pathway=='NABA_MATRISOME_ASSOCIATED',]
top_gene=unlist(top_gene[,8])
NABA_MATR_ASS=as.data.frame(Tcell_1v2DE)
NABA_MATR_ASS=NABA_MATR_ASS[(row.names(NABA_MATR_ASS) %in% top_gene),]
NABA_MATR_ASS=NABA_MATR_ASS[NABA_MATR_ASS$pvalue<0.05,]
if (nrow(NABA_MATR_ASS)!=0){
    NABA_MATR_ASS$pathway="NABA_MATRISOME_ASSOCIATED"
    NABA_MATR_ASS$celltype="Tcell"
    NABA_MATR_ASS$comparison="VEH and TAD"
}
### NABA Tcell 1v2
##  NABA_BASEMENT_MEMBRANE
top_gene=as.data.frame(as.matrix(Tcell_DE$ONEvTWOgs$NABA))
top_gene=top_gene[top_gene$pathway=='NABA_BASEMENT_MEMBRANES',]
top_gene=unlist(top_gene[,8])
NABA_BASE_MEM=as.data.frame(Tcell_1v2DE)
NABA_BASE_MEM=NABA_BASE_MEM[(row.names(NABA_BASE_MEM) %in% top_gene),]
NABA_BASE_MEM=NABA_BASE_MEM[NABA_BASE_MEM$pvalue<0.05,]
if (nrow(NABA_BASE_MEM)!=0){
    NABA_BASE_MEM$pathway="NABA_BASEMENT_MEMBRANE"
    NABA_BASE_MEM$celltype="Tcell"
    NABA_BASE_MEM$comparison="VEH and TAD"
}
### NABA Tcell 1v2
##  NABA_MATRISOME
top_gene=as.data.frame(as.matrix(Tcell_DE$ONEvTWOgs$NABA))
top_gene=top_gene[top_gene$pathway=='NABA_MATRISOME',]
top_gene=unlist(top_gene[,8])
NABA_MAT=as.data.frame(Tcell_1v2DE)
NABA_MAT=NABA_MAT[(row.names(NABA_MAT) %in% top_gene),]
NABA_MAT=NABA_MAT[NABA_MAT$pvalue<0.05,]
if (nrow(NABA_MAT)!=0){
    NABA_MAT$pathway="NABA_MATRISOME"
    NABA_MAT$celltype="Tcell"
    NABA_MAT$comparison="VEH and TAD"
}
### NABA Tcell 1v2
##  NABA_ECM_GLYCOPROTEINS
top_gene=as.data.frame(as.matrix(Tcell_DE$ONEvTWOgs$NABA))
top_gene=top_gene[top_gene$pathway=='NABA_ECM_GLYCOPROTEINS',]
top_gene=unlist(top_gene[,8])
NABA_ECM_GLY=as.data.frame(Tcell_1v2DE)
NABA_ECM_GLY=NABA_ECM_GLY[(row.names(NABA_ECM_GLY) %in% top_gene),]
NABA_ECM_GLY=NABA_ECM_GLY[NABA_ECM_GLY$pvalue<0.05,]
if (nrow(NABA_ECM_GLY)!=0){
    NABA_ECM_GLY$pathway="NABA_ECM_GLYCOPROTEINS"
    NABA_ECM_GLY$celltype="Tcell"
    NABA_ECM_GLY$comparison="VEH and TAD"
}
Tcell1v2 <- rbind(HALL_INFLA_RES,
                  GOMF_EXTRA_tcell,
                  NABA_ECM_REG,
                  NABA_SECRET_FAC,
                  NABA_ECM_AFF,
                  NABA_COLL,
                  NABA_CORE_MAT,
                  NABA_MATR_ASS,
                  NABA_BASE_MEM,
                  NABA_MAT,
                  NABA_ECM_GLY
                 )

meanExpression <- c()
sub = subset(x=sobj,subset=cell_type=="Tcell")
subset = subset(x=sub,subset=group_id==c("VEH","TAD"))
for (i in 1:nrow(Tcell1v2)){
    x = AverageExpression(subset,feature=rownames(Tcell1v2)[i])
    if(length(x)==0){
        meanExpression = c(meanExpression,"NA")
    }
    meanExpression = c(meanExpression,x)
}

g = data.frame(meanExpression = unlist(meanExpression))
Tcell1v2 = cbind(Tcell1v2,g)
write.csv(Tcell1v2,"Tcell_1v2_topGenes.csv")

##### 4v5 Tcell

Tcell4v5 <- data.frame(baseMean=NULL,log2FoldChange=NULL,lfcSE=NULL,stat=NULL,pvalue=NULL,padj=NULL,pathway=NULL,celltype=NULL,comparison=NULL)

### KEGG Tcell 4V5
##  KEGG_MAPK_SIGNALING_PATHWAY
top_gene=as.data.frame(as.matrix(Tcell_DE$FOURvFIVEgs$KEGG))
top_gene=top_gene[top_gene$pathway=='KEGG_MAPK_SIGNALING_PATHWAY',]
top_gene=unlist(top_gene[,8])
KEGG_MAPK=as.data.frame(Tcell_4v5DE)
KEGG_MAPK=KEGG_MAPK[(row.names(KEGG_MAPK) %in% top_gene),]
KEGG_MAPK=KEGG_MAPK[KEGG_MAPK$pvalue<0.05,]
if (nrow(KEGG_MAPK)!=0){
    KEGG_MAPK$pathway="KEGG_MAPK_SIGNALING_PATHWAY"
    KEGG_MAPK$celltype="Tcell"
    KEGG_MAPK$comparison="MLIST and MLIST.TAD"
}
### KEGG Tcell 4V5
##  KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY
top_gene=as.data.frame(as.matrix(Tcell_DE$FOURvFIVEgs$KEGG))
top_gene=top_gene[top_gene$pathway=='KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY',]
top_gene=unlist(top_gene[,8])
KEGG_NOD_LIKE=as.data.frame(Tcell_4v5DE)
KEGG_NOD_LIKE=KEGG_NOD_LIKE[(row.names(KEGG_NOD_LIKE) %in% top_gene),]
KEGG_NOD_LIKE=KEGG_NOD_LIKE[KEGG_NOD_LIKE$pvalue<0.05,]
if (nrow(KEGG_NOD_LIKE)!=0){
    KEGG_NOD_LIKE$pathway="KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY"
    KEGG_NOD_LIKE$celltype="Tcell"
    KEGG_NOD_LIKE$comparison="MLIST and MLIST.TAD"
}
### REACTOME Tcell 4V5
##  REACTOME_INTERLEUKIN_27_SIGNALING
top_gene=as.data.frame(as.matrix(Tcell_DE$FOURvFIVEgs$REACTOME))
top_gene=top_gene[top_gene$pathway=='REACTOME_INTERLEUKIN_27_SIGNALING',]
top_gene=unlist(top_gene[,8])
REACTOME_INTER_27=as.data.frame(Tcell_4v5DE)
REACTOME_INTER_27=REACTOME_INTER_27[(row.names(REACTOME_INTER_27) %in% top_gene),]
REACTOME_INTER_27=REACTOME_INTER_27[REACTOME_INTER_27$pvalue<0.05,]
if (nrow(REACTOME_INTER_27)!=0){
    REACTOME_INTER_27$pathway="REACTOME_INTERLEUKIN_27_SIGNALING"
    REACTOME_INTER_27$celltype="Tcell"
    REACTOME_INTER_27$comparison="MLIST and MLIST.TAD"
}
Tcell4v5 <- rbind(KEGG_MAPK,
                  KEGG_NOD_LIKE,
                  REACTOME_INTER_27)

meanExpression <- c()
sub = subset(x=sobj,subset=cell_type=="Tcell")
subset = subset(x=sub,subset=group_id==c("MLIST","MLIST.TAD"))
for (i in 1:nrow(Tcell4v5)){
    x = AverageExpression(subset,feature=rownames(Tcell4v5)[i])
    if(length(x)==0){
        meanExpression = c(meanExpression,"NA")
    }
    meanExpression = c(meanExpression,x)
}

g = data.frame(meanExpression = unlist(meanExpression))
Tcell4v5 = cbind(Tcell4v5,g)
write.csv(Tcell4v5,"Tcell_4v5_topGenes.csv")
