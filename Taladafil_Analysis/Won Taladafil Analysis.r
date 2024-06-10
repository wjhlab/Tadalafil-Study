library(Seurat)
library(RColorBrewer)
library(SeuratData)
library(dplyr)
library(Seurat)
library(patchwork)
library(MAST)
library(org.Hs.eg.db)

sessionInfo()

sobj <- readRDS("sobj_filtered_transformed_clustered.RDS") # load preprocessed data 

###looking at percent.mt AFTER scaled and mt regressed 
sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")

FeaturePlot(sobj,features = 'percent.mt', label = TRUE, repel = TRUE)

###looking at percent.mt BEFORE scaled and mt regressed 
FeaturePlot(sobj,features = 'percent.mt', label = TRUE, repel = TRUE)

levels(sobj)

new.cluster.ids <- c("Myeloid", "Granulocyte", "Tumor", "Macrophage", "Tcell", "CAF",
    "DC", "Endothelial")

names(new.cluster.ids) <- levels(sobj)
#subseurat <- RenameIdents(subseurat, new.cluster.ids)
#DimPlot(subseurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

levels(sobj)

pdf('sample_umap.pdf',width=15,height=15)
DimPlot(subseurat, reduction = "umap",label=FALSE,group.by = 'sample_id')
dev.off()

sobj@assays

DimPlot(subseurat, reduction = "umap",label=FALSE,group.by = 'sample_id')

#Use FindAllMarkers to get top DE genes between clusters for annotation

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
sobj.markers <- FindAllMarkers(sobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

sobj.markers<- ssobj.markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)

unique(levels(sobj))

#B_cells<- c("Cd19","Cd79b","Cd45r","Cd79a","B220")
#T_cells<- c("Cd3d","Cd3e","Cd3g")
#Neutrophils<-c("Itgam","S100a9","Cxcr4","Ly6g")
#Macrophage<- c("Adgre1", "Cd68", "Ptprc", "Fcgr1", "Mertk", "Itgam", "Csf1r", "Ccr2")
#M1<- c("Ccl2","Cxcl1","Cd163","Wfdc17","Tnfa","Tnfsf1a","Tnf")
#M2<- c("Arg1","Ccl22","Ccl18","Stat6","Cd206", "Mrc1")


#Myeloid
#FeaturePlot(sobj, features=c("Arg1","C1qc","H2-Aa","H2-Eb1","Lgmn","Ctsl"))

#T cell
#FeaturePlot(sobj, features = c("Cd3d","Cd3e","Cd3g"))

#Granulocyte
#FeaturePlot(sobj,features=c("Itgam","S100a9","Cxcr4"))

#Fibroblasts
#FeaturePlot(sobj,features=c("Dcn", "Col12a1", "Mmp2"))

#Macs
#FeaturePlot(sobj,features=c("Adgre1", "Cd68", "Fcgr1", "Itgam", "Csf1r", "Ccr2","Mrc1"))

#M1
#FeaturePlot(sobj,features=c("Ccl2","Cxcl1","Cd163","Wfdc17","Tnfa","Tnfsf1a","Tnf"))

#M2
#FeaturePlot(sobj,features=c("Arg1","Ccl22","Ccl18","Stat6", "Mrc1"))

#FeaturePlot(sobj,features=c("Cyp11a1","Il4","Itga2b","Cd200r3","Slc18a2","Il1rl1","Ccl4","Gzmb"))

#FeaturePlot(sobj,features=c("S100a8", "S100a9", "S100a12"))

#B cells
#FeaturePlot(sobj,features=c("Igkc", "Cd55", "Fcmr", "Ms4a1"))

#DCs
#FeaturePlot(sobj,features=c("Flt3", "Batf3"))

#NK
#FeaturePlot(sobj,features=c("Nkg7","Klrb1", "Klrk1", "Ncr1", "Klrc1", "Cd56", "Fcgr3a", "Fcg3","Fcer1g"))

#monocyte
#FeaturePlot(sobj,features=c("Cd14","Itgal","Csf1r","Cd44","Ccl2","Ccr2","Itgam","Fcgr3","Ly6C1"))

#Endothelial
FeaturePlot(sobj,features=c("Pecam1","Vwf"))

#Cancer/epithelial
#FeaturePlot(sobj,features=c("Wfdc2", "Krt8", "Spp1","Krt18"))

new.cluster.ids <- c("Myeloid", "Granulocyte", "Tumor", "Macrophage", "T_cell", "CAF",
    "DC", "Endothelial")

names(new.cluster.ids) <- levels(sobj)
sobj <- RenameIdents(sobj, new.cluster.ids)
DimPlot(sobj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(sobj, file = "Annotated_clusters.RDS")

#Load annotated data
sobj <- readRDS("Annotated_clusters.RDS") 

head(sobj@meta.data)

#Subset seurat into different cell type objects
sobj_gra = subset(x =sobj, subset = seurat_clusters == 1)
sobj_Tcell = subset(x =sobj, subset = seurat_clusters == 4)
sobj_macro = subset(x =sobj, subset = seurat_clusters == 3)
sobj_DC = subset(x =sobj, subset = seurat_clusters == 6)

#Look at cell cycle states of cells

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

library(org.Hs.eg.db)
library(homologene)
# Convert human gene symbol to mouse ortholog symbol since function uses human genes
convert_human_to_mouse <- function(human_gene) {
  mouse_gene <- human2mouse(human_gene)
  return(mouse_gene)
}
# Example usage
human_gene <- s.genes
mouse_ortholog <- convert_human_to_mouse(human_gene)
print(mouse_ortholog)

s.genes.mouse<-mouse_ortholog

human_gene <- g2m.genes
mouse_ortholog2 <- convert_human_to_mouse(human_gene)
print(mouse_ortholog)

g2m.genes.mouse<-mouse_ortholog2

mouse_genes <- rownames(sobj)

s.genes.mouse$mouseGene

sobj@meta.data

#score cells based on cell cycle state
sobj <- CellCycleScoring(sobj, s.features = s.genes.mouse$mouseGene, g2m.features = g2m.genes.mouse$mouseGene, set.ident = TRUE)

head(sobj@meta.data)

#plot UMAP with cell cycle scores 
sobj <- RunPCA(sobj, features = c(s.genes.mouse$mouseGene, g2m.genes.mouse$mouseGene))
DimPlot(sobj)

FeaturePlot(sobj, features = c("percent.mt","nCount_RNA","G2M.Score","S.Score"))



sample_id = list("G1mus1","G1mus2","G1mus3","G2mus1","G2mus2","G2mus3","G2mus3_rep","G3mus1","G3mus2","G3mus3","G4mus1","G4mus2","G4mus3","G5mus1","G5mus2","G5mus3")
tube_id = list('1.1','1.2','1.3','2.1','2.2','2.3','2.3.1','3.1','3.2','3.3','4.1','4.2','4.3','5.1','5.2','5.3')
group = list('VEH','VEH','VEH','TAD','TAD','TAD','TAD','CLIST','CLIST','CLIST','MLIST','MLIST','MLIST','MLIST_TAD','MLIST_TAD','MLIST_TAD')
index_id = list('B1','B2','B3','B4','B5','B6','C4','B7','B8','B9','B10','B11','B12','C1','C2','C3')

head(sobj@meta.data)

sobj@meta.data['group_id']=NA
sobj@meta.data['sample_id']=NA
sobj@meta.data['tube_id']=NA
sobj@meta.data['index']=NA

head(sobj@meta.data)

for (i in 1:nrow(sobj@meta.data)){
    indicator = rownames(sobj@meta.data)[i]
    indicator = gsub(".*-","",indicator)
    sobj@meta.data[i,10]=group[as.numeric(indicator)]
    sobj@meta.data[i,11]=sample_id[as.numeric(indicator)]
    sobj@meta.data[i,12]=tube_id[as.numeric(indicator)]
    sobj@meta.data[i,13]=index_id[as.numeric(indicator)]
}

unique(levels(sobj))



sobj_backup = sobj

sobj_backup@meta.data





#Pathway gene sets
load("hallmark_pathways.rda")
load("kegg_pathways.rda")
HALL_INFLA_RES=hallmark_pathways$HALLMARK_INFLAMMATORY_RESPONSE
HALL_INFLA_RES <- convert_human_to_mouse(HALL_INFLA_RES)$mouseGene

HALL_ANGI=hallmark_pathways$HALLMARK_ANGIOGENESIS
HALL_ANGI <- convert_human_to_mouse(HALL_ANGI)$mouseGene

HALL_REA_O2=hallmark_pathways$HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY
HALL_REA_O2 <- convert_human_to_mouse(HALL_REA_O2)$mouseGene

HALL_APOP=hallmark_pathways$HALLMARK_APOPTOSIS
HALL_APOP <- convert_human_to_mouse(HALL_APOP)$mouseGene

mmu00330=kegg$`mmu00330 Arginine and proline metabolism`
mmu00330 <- convert_human_to_mouse(mmu00330)$mouseGene

mmu04022=kegg$`mmu04022 cGMP-PKG signaling pathway`
mmu04022 <- convert_human_to_mouse(mmu04022)$mouseGene

mmu04060=kegg$`mmu04060 Cytokine-cytokine receptor interaction`
mmu04060 <- convert_human_to_mouse(mmu04060)$mouseGene

mmu04062=kegg$`mmu04062 Chemokine signaling pathway`
mmu04062 <- convert_human_to_mouse(mmu04062)$mouseGene

mmu04150=kegg$`mmu04150 mTOR signaling pathway`
mmu04150 <- convert_human_to_mouse(mmu04150)$mouseGene

mmu04218=kegg$`mmu04218 Cellular senescence`
mmu04218 <- convert_human_to_mouse(mmu04218)$mouseGene

mmu04370=kegg$`mmu04370 VEGF signaling pathway`
mmu04370 <- convert_human_to_mouse(mmu04370)$mouseGene

mmu04512=kegg$`mmu04512 ECM-receptor interaction`
mmu04512 <- convert_human_to_mouse(mmu04512)$mouseGene

mmu04612=kegg$`mmu04612 Antigen processing and presentation`
mmu04612 <- convert_human_to_mouse(mmu04612)$mouseGene

mmu05235=kegg$`mmu05235 PD-L1 expression and PD-1 checkpoint pathway in cancer`
mmu05235 <- convert_human_to_mouse(mmu05235)$mouseGene


geneSetsList1 <- list(Tcells = c("Cd2","Cd3d","Cd3e","Cd3g","Cd28","Cd4","Cd8b1","Foxp3","Itgae","Cd3e", "Cd4","Foxp3","Itgae","Cd25","Il2r","Il-2r", "Il2ra", "Il-2ra","Hla","Cd38","Cd39","Cd69","Tcra","Tcrb","Tcr","H2-eb1","H2eb1","Lag3","Cd244","Eomes","Ptger4","Cd44","Cd62l","Cd45ra", "Ccr7", "Cd62l", "Cd127", "Cd132"),
                        Bcells = c("Cd19","Cd79b","Cd45r","Cd79a","B220","CD19","B220", "CD45R", "CD45r"),
                        Macro = c("F4/80","Adgre1","Cd68","Cd45","Ptprc","Cd64","Fcgr1","Mertk","Ccl2","Cxcl1","Cd163","Wfdc17","Tnfa","Tnfsf1a","Tnf","Arg1","Ccl22","Ccl18","Stat6","Cd206", "Mrc1"),
                        neutrophil_genes <-c("Itgam","S100a9","Cxcr4","Ly6g"),
                        monocyte_macrophage_genes <- c("Adgre1", "Cd68", "Ptprc", "Fcgr1", "Mertk", "Itgam", "Csf1r", "Ccr2"),
                        M1 <- c("Ccl2","Cxcl1","Cd163","Wfdc17","Tnfa","Tnfsf1a","Tnf"),
                        M2 <- c("Arg1","Ccl22","Ccl18","Stat6","Cd206", "Mrc1"))

geneSetsList2 <- list(HALL_INFLA_RES,
                      HALL_ANGI,
                      HALL_REA_O2,
                      HALL_APOP,
                      mmu00330,
                      mmu04022,
                      mmu04060,
                      mmu04062,
                      mmu04150,
                      mmu04218,
                      mmu04370,
                      mmu04512,
                      mmu04612,
                      mmu05235
)



sobj <- AddModuleScore(sobj, geneSetsList2, pool = NULL, nbin = 24, ctrl = 100,
                       k = FALSE, assay = NULL, name = "module", seed = 1)

colnames(sobj@meta.data)

names(sobj@meta.data)[names(sobj@meta.data) == 'module1'] <- 'HALL_INFLA_RES'
names(sobj@meta.data)[names(sobj@meta.data) == 'module2'] <- 'HALL_ANGI'
names(sobj@meta.data)[names(sobj@meta.data) == 'module3'] <- 'HALL_REA_O2'
names(sobj@meta.data)[names(sobj@meta.data) == 'module4'] <- 'HALL_APOP'
names(sobj@meta.data)[names(sobj@meta.data) == 'module5'] <- 'mmu00330'
names(sobj@meta.data)[names(sobj@meta.data) == 'module6'] <- 'mmu04022'
names(sobj@meta.data)[names(sobj@meta.data) == 'module7'] <- 'mmu04060'
names(sobj@meta.data)[names(sobj@meta.data) == 'module8'] <- 'mmu04062'
names(sobj@meta.data)[names(sobj@meta.data) == 'module9'] <- 'mmu04150'
names(sobj@meta.data)[names(sobj@meta.data) == 'module10'] <- 'mmu04218'
names(sobj@meta.data)[names(sobj@meta.data) == 'module11'] <- 'mmu04370'
names(sobj@meta.data)[names(sobj@meta.data) == 'module12'] <- 'mmu04512'
names(sobj@meta.data)[names(sobj@meta.data) == 'module13'] <- 'mmu04612'
names(sobj@meta.data)[names(sobj@meta.data) == 'module14'] <- 'mmu05235'

levels(sobj)



##### Carter's Codes for Violin plots 2/10/2023

sobj <- readRDS("Annotated_clusters.RDS") 

library(org.Hs.eg.db)
library(homologene)

#Look at cell cycle states of cells

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Convert human gene symbol to mouse ortholog symbol since function uses human genes
convert_human_to_mouse <- function(human_gene) {
  mouse_gene <- human2mouse(human_gene)
  return(mouse_gene)
}
# Example usage
human_gene <- s.genes
mouse_ortholog <- convert_human_to_mouse(human_gene)

s.genes.mouse<-mouse_ortholog

human_gene <- g2m.genes
mouse_ortholog2 <- convert_human_to_mouse(human_gene)
g2m.genes.mouse<-mouse_ortholog2

sample_id = list("G1mus1","G1mus2","G1mus3","G2mus1","G2mus2","G2mus3","G2mus3_rep","G3mus1","G3mus2","G3mus3","G4mus1","G4mus2","G4mus3","G5mus1","G5mus2","G5mus3")
tube_id = list('1.1','1.2','1.3','2.1','2.2','2.3','2.3.1','3.1','3.2','3.3','4.1','4.2','4.3','5.1','5.2','5.3')
group = list('VEH','VEH','VEH','TAD','TAD','TAD','TAD','CLIST','CLIST','CLIST','MLIST','MLIST','MLIST','MLIST_TAD','MLIST_TAD','MLIST_TAD')
index_id = list('B1','B2','B3','B4','B5','B6','C4','B7','B8','B9','B10','B11','B12','C1','C2','C3')

sobj@meta.data['group_id']=NA
sobj@meta.data['sample_id']=NA
sobj@meta.data['tube_id']=NA
sobj@meta.data['index']=NA

head(sobj@meta.data)



for (i in 1:nrow(sobj@meta.data)){
    indicator = rownames(sobj@meta.data)[i]
    indicator = gsub(".*-","",indicator)
    sobj@meta.data[i,10]=group[as.numeric(indicator)]
    sobj@meta.data[i,11]=sample_id[as.numeric(indicator)]
    sobj@meta.data[i,12]=tube_id[as.numeric(indicator)]
    sobj@meta.data[i,13]=index_id[as.numeric(indicator)]
}

#DimPlot(sobj, reduction = "umap",label=FALSE,group.by = 'group_id')

sobj_backup=sobj

#Pathway gene sets
load("hallmark_pathways.rda")
load("kegg_pathways.rda")
HALL_INFLA_RES=hallmark_pathways$HALLMARK_INFLAMMATORY_RESPONSE
HALL_INFLA_RES <- convert_human_to_mouse(HALL_INFLA_RES)$mouseGene

HALL_ANGI=hallmark_pathways$HALLMARK_ANGIOGENESIS
HALL_ANGI <- convert_human_to_mouse(HALL_ANGI)$mouseGene

HALL_REA_O2=hallmark_pathways$HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY
HALL_REA_O2 <- convert_human_to_mouse(HALL_REA_O2)$mouseGene

HALL_APOP=hallmark_pathways$HALLMARK_APOPTOSIS
HALL_APOP <- convert_human_to_mouse(HALL_APOP)$mouseGene

Arginine_proline=kegg$`mmu00330 Arginine and proline metabolism`
Arginine_proline <- convert_human_to_mouse(Arginine_proline)$mouseGene

cGMP_PKG=kegg$`mmu04022 cGMP-PKG signaling pathway`
cGMP_PKG <- convert_human_to_mouse(cGMP_PKG)$mouseGene

Cytokine_cytokine=kegg$`mmu04060 Cytokine-cytokine receptor interaction`
Cytokine_cytokine <- convert_human_to_mouse(Cytokine_cytokine)$mouseGene

Chemokine=kegg$`mmu04062 Chemokine signaling pathway`
Chemokine <- convert_human_to_mouse(Chemokine)$mouseGene

mTOR=kegg$`mmu04150 mTOR signaling pathway`
mTOR <- convert_human_to_mouse(mTOR)$mouseGene

Cellular_Senescence=kegg$`mmu04218 Cellular senescence`
Cellular_Senescence <- convert_human_to_mouse(Cellular_Senescence)$mouseGene

VEGF=kegg$`mmu04370 VEGF signaling pathway`
VEGF <- convert_human_to_mouse(VEGF)$mouseGene

ECM_receptor=kegg$`mmu04512 ECM-receptor interaction`
ECM_receptor <- convert_human_to_mouse(ECM_receptor)$mouseGene

Antigen=kegg$`mmu04612 Antigen processing and presentation`
Antigen <- convert_human_to_mouse(Antigen)$mouseGene

PD_L1_PD_1=kegg$`mmu05235 PD-L1 expression and PD-1 checkpoint pathway in cancer`
PD_L1_PD_1 <- convert_human_to_mouse(PD_L1_PD_1)$mouseGene


geneSetsList1 <- list(Tcells = c("Cd2","Cd3d","Cd3e","Cd3g","Cd28","Cd4","Cd8b1","Foxp3","Itgae","Cd3e", "Cd4","Foxp3","Itgae","Cd25","Il2r","Il-2r", "Il2ra", "Il-2ra","Hla","Cd38","Cd39","Cd69","Tcra","Tcrb","Tcr","H2-eb1","H2eb1","Lag3","Cd244","Eomes","Ptger4","Cd44","Cd62l","Cd45ra", "Ccr7", "Cd62l", "Cd127", "Cd132"),
                        Bcells = c("Cd19","Cd79b","Cd45r","Cd79a","B220","CD19","B220", "CD45R", "CD45r"),
                        Macro = c("F4/80","Adgre1","Cd68","Cd45","Ptprc","Cd64","Fcgr1","Mertk","Ccl2","Cxcl1","Cd163","Wfdc17","Tnfa","Tnfsf1a","Tnf","Arg1","Ccl22","Ccl18","Stat6","Cd206", "Mrc1"),
                        neutrophil_genes <-c("Itgam","S100a9","Cxcr4","Ly6g"),
                        monocyte_macrophage_genes <- c("Adgre1", "Cd68", "Ptprc", "Fcgr1", "Mertk", "Itgam", "Csf1r", "Ccr2"),
                        M1 <- c("Ccl2","Cxcl1","Cd163","Wfdc17","Tnfa","Tnfsf1a","Tnf"),
                        M2 <- c("Arg1","Ccl22","Ccl18","Stat6","Cd206", "Mrc1"))

geneSetsList2 <- list(HALL_INFLA_RES,
                      HALL_ANGI,
                      HALL_REA_O2, 
                      HALL_APOP,
                      Arginine_proline,
                      cGMP_PKG,
                      Cytokine_cytokine,
                      Chemokine,
                      mTOR,
                      Cellular_Senescence,
                      VEGF,
                      ECM_receptor,
                      Antigen,
                      PD_L1_PD_1
)



sobj <- AddModuleScore(sobj, geneSetsList2, pool = NULL, nbin = 24, ctrl = 100,
                       k = FALSE, assay = NULL, name = "module", seed = 1)

names(sobj@meta.data)[names(sobj@meta.data) == 'module1'] <- 'HALL_INFLA_RES'
names(sobj@meta.data)[names(sobj@meta.data) == 'module2'] <- 'HALL_ANGI'
names(sobj@meta.data)[names(sobj@meta.data) == 'module3'] <- 'HALL_REA_O2'
names(sobj@meta.data)[names(sobj@meta.data) == 'module4'] <- 'HALL_APOP'
names(sobj@meta.data)[names(sobj@meta.data) == 'module5'] <- 'Arginine_proline'
names(sobj@meta.data)[names(sobj@meta.data) == 'module6'] <- 'cGMP_PKG'
names(sobj@meta.data)[names(sobj@meta.data) == 'module7'] <- 'Cytokine_cytokine'
names(sobj@meta.data)[names(sobj@meta.data) == 'module8'] <- 'Chemokine'
names(sobj@meta.data)[names(sobj@meta.data) == 'module9'] <- 'mTOR'
names(sobj@meta.data)[names(sobj@meta.data) == 'module10'] <- 'Cellular senescence'
names(sobj@meta.data)[names(sobj@meta.data) == 'module11'] <- 'VEGF'
names(sobj@meta.data)[names(sobj@meta.data) == 'module12'] <- 'ECM_receptor'
names(sobj@meta.data)[names(sobj@meta.data) == 'module13'] <- 'Antigen'
names(sobj@meta.data)[names(sobj@meta.data) == 'module14'] <- 'PD_L1_PD_1'

colnames(sobj@meta.data)



library(ggplot2)

sub_cluster0 = subset(x =sobj, subset = seurat_clusters == 0)
plot0=VlnPlot(sub_cluster0, features = 'Antigen', group.by = "group_id",slot="scale.data",pt.size = 0)+labs(subtitle='Myeloid')+stat_summary(fun.y = median, geom='crossbar', size = .3, colour = "black")

sub_cluster1 = subset(x =sobj, subset = seurat_clusters == 1)
plot1=VlnPlot(sub_cluster1, features = 'Antigen', group.by = "group_id",slot="scale.data",pt.size = 0)+labs(subtitle='Granulocytes')+stat_summary(fun.y = median, geom='crossbar', size = .3, colour = "black")

sub_cluster2 = subset(x =sobj, subset = seurat_clusters == 2)
plot2=VlnPlot(sub_cluster2, features = 'Antigen', group.by = "group_id",slot="scale.data",pt.size = 0)+labs(subtitle='Tumor')+stat_summary(fun.y = median, geom='crossbar', size = .3, colour = "black")

sub_cluster3 = subset(x =sobj, subset = seurat_clusters == 3)
plot3=VlnPlot(sub_cluster3, features = 'Antigen', group.by = "group_id",slot="scale.data",pt.size = 0)+labs(subtitle='Macrophage')+stat_summary(fun.y = median, geom='crossbar', size = .3, colour = "black")

sub_cluster4 = subset(x =sobj, subset = seurat_clusters == 4)
plot4=VlnPlot(sub_cluster4, features = 'Antigen', group.by = "group_id",slot="scale.data",pt.size = 0)+labs(subtitle='T cell')+stat_summary(fun.y = median, geom='crossbar', size = .3, colour = "black")

sub_cluster5 = subset(x =sobj, subset = seurat_clusters == 5)
plot5=VlnPlot(sub_cluster5, features = 'Antigen', group.by = "group_id",slot="scale.data",pt.size = 0)+labs(subtitle='CAF')+stat_summary(fun.y = median, geom='crossbar', size = .3, colour = "black")

sub_cluster6 = subset(x =sobj, subset = seurat_clusters == 6)
plot6=VlnPlot(sub_cluster6, features = 'Antigen', group.by = "group_id",slot="scale.data",pt.size = 0)+labs(subtitle='DC')+stat_summary(fun.y = median, geom='crossbar', size = .3, colour = "black")

sub_cluster7 = subset(x =sobj, subset = seurat_clusters == 7)
plot7=VlnPlot(sub_cluster7, features = 'Antigen', group.by = "group_id",slot="scale.data",pt.size = 0)+labs(subtitle='Monocyte')+stat_summary(fun.y = median, geom='crossbar', size = .3, colour = "black")

pdf('9.pdf',width=20,height=10)
CombinePlots(ncol=4,
plots = list(plot0, plot1, plot2, plot3, plot4, plot5, plot6, plot7)
)
dev.off()


