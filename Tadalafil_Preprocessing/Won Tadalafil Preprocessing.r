library(Seurat)
library(RColorBrewer)
library(SeuratData)
library(dplyr)
library(Seurat)
library(patchwork)
library(MAST)

# load unfiltered data 
sobj <- readRDS("unfiltered_sobj.RDS") 

file.info("unfiltered_sobj.RDS")

#saveRDS(sobj, file ="unfiltered_sobj.RDS")

# store mitochondrial percentage in object meta data
sobj <- PercentageFeatureSet(sobj, pattern = "^mt-", col.name = "percent.mt")

# run sctransform
# sobj <- SCTransform(sobj, vars.to.regress = "percent.mt", verbose = FALSE)

# qc plots
VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# filtering based on above plots
sobj <- subset(sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500  & percent.mt < 5)

#Normalize data
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)

#ID highly variable genes (feature selection)
sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sobj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sobj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scale based on mito gene expression 
#to prevent UMAP from being driven by mt gene expression
#this takes 4 hours to run
all.genes <- rownames(sobj)
sobj <- ScaleData(object = sobj,features=all.genes, vars.to.regress = c("percent.mt"))


saveRDS(sobj, file = "sobj_filtered_transformed_scaled.RDS")

# load data 
sobj <- readRDS("sobj_filtered_transformed_scaled.RDS") 

#Linear dimensional reduction
sobj <- RunPCA(sobj, features = VariableFeatures(object = sobj))

# Examine and visualize PCA results a few different ways
print(sobj[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(sobj, dims = 1:2, reduction = "pca")

DimPlot(sobj, reduction = "pca")

DimHeatmap(sobj, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(sobj, dims = 1:15, cells = 500, balanced = TRUE)

colnames(sobj)

#determine dimensionality of data
ElbowPlot(sobj)

#Cluster cells
#dimensions determine how many/where dots in the UMAP occurs
#resolution determines how many groups/clusters (colors) are in the UMAP
sobj <- FindNeighbors(sobj, dims = 1:50)
sobj <- FindClusters(sobj, resolution = 0.2)

#reticulate::py_install(packages ='umap-learn')

#n.neighbors determines how many cells are needed to make a neighborhood
sobj <- RunUMAP(sobj, dims = 1:50, n.neighbors = 50L)

#Plot the UMAP
DimPlot(sobj, reduction = "umap")

#clean up the UMAP
remove <- c(rownames(sobj@meta.data[sobj$seurat_clusters==16,]), #remove clusters that have less than ~150 cells
            rownames(sobj@meta.data[sobj$seurat_clusters=="singleton",])) #remove singletons
            

#Making object that has small clusters & singletons removed
keep <- setdiff(colnames(sobj),remove)
length(keep)
dim(sobj)

sobj <- sobj[,keep]

sobj <- FindNeighbors(sobj, dims = 1:50)
sobj <- FindClusters(sobj, resolution = 0.05)
DimPlot(sobj, reduction = "umap")

#Look at mito gene expression across clusters
FeaturePlot(sobj,features = 'percent.mt', label = TRUE, repel = TRUE)

DimPlot(sobj, reduction = "umap")

saveRDS(sobj, file = "sobj_filtered_transformed_clustered.RDS")
