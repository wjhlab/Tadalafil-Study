library(domino)
library(Seurat)
library(SeuratObject)
library(SeuratObject)

DominoInteractions <- setClass(
  Class = "DominoInteractions",
  slots = c(
    subjects = "factor", # identification for the different domino results
    subject_meta = "data.frame", # decoding data frame for covariates
    subject_features = "list" # nested lists of:
    # subject name
    # cluster
    # cluster tfs
    # cluster receptors
    # incoming ligands
  )
)


re_list <- readRDS(file = "E:/WJHLab/revised_domino/tadalafil_domino_features.rds")
dom_clist <- readRDS(file = "E:/WJHLab/revised_domino/new_unbuilt/J1568_CLIST_domino_unbuilt.rds")
dom_mlist <- readRDS(file = "E:/WJHLab/revised_domino/new_unbuilt/J1568_MLIST_domino_unbuilt.rds")
dom_mlisttad <- readRDS(file = "E:/WJHLab/revised_domino/new_unbuilt/J1568_MLISTTAD_domino_unbuilt.rds")
dom_tad <- readRDS(file = "E:/WJHLab/revised_domino/new_unbuilt/J1568_TAD_domino_unbuilt.rds")
dom_veh <- readRDS(file = "E:/WJHLab/revised_domino/new_unbuilt/J1568_VEH_domino_unbuilt.rds")


count_interactions <- 
  function(DominoInteractions, cluster, interaction = "rec_lig", group.by = NULL,
           subject_id = "PubID"){
    di <- DominoInteractions
    
    subjects <- di@subjects
    all_int <- c()
    for(id in subjects){
      interact <- di@subject_features[[id]][[cluster]][[interaction]]
      all_int <- c(all_int, interact)
    }
    df <- as.data.frame(table(all_int))
    rownames(df) <- df[["all_int"]]
    
    if(!is.null(group.by)){
      if(!group.by %in% colnames(di@subject_meta)){
        stop("group.by variable not present in subject_meta")
      }
      groups <- levels(factor(di@subject_meta[[group.by]]))
      
      for(g in groups){
        g_index <- di@subject_meta[[group.by]] == g
        g_subjects <- di@subject_meta[g_index, subject_id]
        
        # count occurrence of these interactions with grep
        int_count <- list()
        for(i in rownames(df)){
          count<- sapply(g_subjects, function(x){
            i %in% di@subject_features[[x]][[cluster]][[interaction]]
          })
          int_count[[i]] <-
            sum(count)
        }
        df[[g]] <- int_count
      }
    }
    return(df)
  }


comp_table <- data.frame(all_int = NULL,
                         Freq = NULL,
                         CLIST = NULL,
                         MLIST = NULL,
                         MLIST.TAD = NULL,
                         TAD = NULL,
                         VEH = NULL,
                         cell_type = NULL)
cell_type <- c("CAF","DC","Endothelial","Granulocyte","Macrophage","Myeloid","Tcell","Tumor")

for (i in 1:8){
  comparison <- count_interactions(re_list,
                                   cluster = cell_type[i], 
                                   interaction = "rec",
                                   group.by = "Treatment",
                                   subject_id = "Treatment")
  comparison$cell_type <- c(rep(cell_type[i],nrow(comparison)))
  comp_table <- rbind(comp_table,comparison)
}

############for all five group setting
#setwd("E:/WJHLab/domino/meeting_3")
veh1tad0 <- comp_table[comp_table$VEH==1&comp_table$TAD==0,]
tad1veh0 <- comp_table[comp_table$VEH==0&comp_table$TAD==1,]
mlist1mlisttad0 <- comp_table[comp_table$MLIST==1&comp_table$MLISTTAD==0,]
mlisttad1mlist0 <- comp_table[comp_table$MLIST==0&comp_table$MLISTTAD==1,]
tadmlisttad <- comp_table[comp_table$TAD==1&comp_table$MLISTTAD==1&comp_table$MLIST==0&comp_table$CLIST==0&comp_table$VEH==0,]
non_tad <- comp_table[comp_table$VEH==1&comp_table$CLIST==1&comp_table$MLIST==1,]
df <- data.frame(lapply(non_tad, as.character), stringsAsFactors = FALSE)
#write.csv(df,file="non_tad.csv")

dim(veh1tad0)
dim(tad1veh0)
dim(mlist1mlisttad0)
dim(mlisttad1mlist0)
dim(tadmlisttad)


veh1tad0 <- veh1tad0[,!colnames(veh1tad0)%in%c("Freq")]
colnames(veh1tad0)[which(names(veh1tad0) == "all_int")] <- "receptors"
receptor_tem <- c()
cell_type_tem <- c()
ligand <- c()
for (i in 1:nrow(veh1tad0)){
  rec <- as.character(veh1tad0[i,1])
  cell <- as.character(veh1tad0[i,7])
  temp <- dom_veh@linkages$rec_lig[[rec]]
  ligand <- c(ligand,temp)
  receptor_tem <- c(receptor_tem,rep(rec,length(temp)))
  cell_type_tem <- c(cell_type_tem,rep(cell,length(temp)))
}
new_veh1tad0 <- data.frame(receptor = receptor_tem,
                    CLIST = rep(0,length(receptor_tem)),
                    MLIST = rep(0,length(receptor_tem)),
                    MLIST.TAD = rep(0,length(receptor_tem)),
                    TAD = rep(0,length(receptor_tem)),
                    VEH = rep(1,length(receptor_tem)),
                    cell_type = cell_type_tem,
                    ligand = ligand)
#saveRDS(new_veh1tad0,file = 'veh1tad0.rds')

tad1veh0 <- tad1veh0[,!colnames(tad1veh0)%in%c("Freq")]
colnames(tad1veh0)[which(names(tad1veh0) == "all_int")] <- "receptors"
receptor_tem <- c()
cell_type_tem <- c()
ligand <- c()
for (i in 1:nrow(tad1veh0)){
  rec <- as.character(tad1veh0[i,1])
  cell <- as.character(tad1veh0[i,7])
  temp <- dom_tad@linkages$rec_lig[[rec]]
  ligand <- c(ligand,temp)
  receptor_tem <- c(receptor_tem,rep(rec,length(temp)))
  cell_type_tem <- c(cell_type_tem,rep(cell,length(temp)))
}
new_tad1veh0 <- data.frame(receptor = receptor_tem,
                           CLIST = rep(0,length(receptor_tem)),
                           MLIST = rep(0,length(receptor_tem)),
                           MLIST.TAD = rep(0,length(receptor_tem)),
                           TAD = rep(0,length(receptor_tem)),
                           VEH = rep(1,length(receptor_tem)),
                           cell_type = cell_type_tem,
                           ligand = ligand)
#saveRDS(new_tad1veh0,file = 'tad1veh0.rds')

mlist1mlisttad0 <- mlist1mlisttad0[,!colnames(mlist1mlisttad0)%in%c("Freq")]
colnames(mlist1mlisttad0)[which(names(mlist1mlisttad0) == "all_int")] <- "receptors"
receptor_tem <- c()
cell_type_tem <- c()
ligand <- c()
for (i in 1:nrow(mlist1mlisttad0)){
  rec <- as.character(mlist1mlisttad0[i,1])
  cell <- as.character(mlist1mlisttad0[i,7])
  temp <- dom_tad@linkages$rec_lig[[rec]]
  ligand <- c(ligand,temp)
  receptor_tem <- c(receptor_tem,rep(rec,length(temp)))
  cell_type_tem <- c(cell_type_tem,rep(cell,length(temp)))
}
new_mlist1mlisttad0 <- data.frame(receptor = receptor_tem,
                           CLIST = rep(0,length(receptor_tem)),
                           MLIST = rep(0,length(receptor_tem)),
                           MLIST.TAD = rep(0,length(receptor_tem)),
                           TAD = rep(0,length(receptor_tem)),
                           VEH = rep(1,length(receptor_tem)),
                           cell_type = cell_type_tem,
                           ligand = ligand)
#saveRDS(new_mlist1mlisttad0,file = 'new_mlist1mlisttad0.rds')

mlisttad1mlist0 <- mlisttad1mlist0[,!colnames(mlisttad1mlist0)%in%c("Freq")]
colnames(mlisttad1mlist0)[which(names(mlisttad1mlist0) == "all_int")] <- "receptors"
receptor_tem <- c()
cell_type_tem <- c()
ligand <- c()
for (i in 1:nrow(mlisttad1mlist0)){
  rec <- as.character(mlisttad1mlist0[i,1])
  cell <- as.character(mlisttad1mlist0[i,7])
  temp <- dom_tad@linkages$rec_lig[[rec]]
  ligand <- c(ligand,temp)
  receptor_tem <- c(receptor_tem,rep(rec,length(temp)))
  cell_type_tem <- c(cell_type_tem,rep(cell,length(temp)))
}
new_mlisttad1mlist0 <- data.frame(receptor = receptor_tem,
                                  CLIST = rep(0,length(receptor_tem)),
                                  MLIST = rep(0,length(receptor_tem)),
                                  MLIST.TAD = rep(0,length(receptor_tem)),
                                  TAD = rep(0,length(receptor_tem)),
                                  VEH = rep(1,length(receptor_tem)),
                                  cell_type = cell_type_tem,
                                  ligand = ligand)
#saveRDS(new_mlisttad1mlist0,file = 'new_mlisttad1mlist0.rds')

tadmlisttad <- tadmlisttad[,!colnames(tadmlisttad)%in%c("Freq")]
colnames(tadmlisttad)[which(names(tadmlisttad) == "all_int")] <- "receptors"
receptor_tem <- c()
cell_type_tem <- c()
ligand <- c()
for (i in 1:nrow(tadmlisttad)){
  rec <- as.character(tadmlisttad[i,1])
  cell <- as.character(tadmlisttad[i,7])
  temp <- dom_tad@linkages$rec_lig[[rec]]
  ligand <- c(ligand,temp)
  receptor_tem <- c(receptor_tem,rep(rec,length(temp)))
  cell_type_tem <- c(cell_type_tem,rep(cell,length(temp)))
}
new_tadmlisttad <- data.frame(receptor = receptor_tem,
                                  CLIST = rep(0,length(receptor_tem)),
                                  MLIST = rep(0,length(receptor_tem)),
                                  MLIST.TAD = rep(0,length(receptor_tem)),
                                  TAD = rep(0,length(receptor_tem)),
                                  VEH = rep(1,length(receptor_tem)),
                                  cell_type = cell_type_tem,
                                  ligand = ligand)
#saveRDS(new_tadmlisttad,file = 'tadmlisttad.RDS')

non_tad <- non_tad[,!colnames(non_tad)%in%c("Freq")]
colnames(non_tad)[which(names(non_tad) == "all_int")] <- "receptors"
receptor_tem <- c()
cell_type_tem <- c()
ligand <- c()
for (i in 1:nrow(non_tad)){
  rec <- as.character(non_tad[i,1])
  cell <- as.character(non_tad[i,7])
  temp <- dom_tad@linkages$rec_lig[[rec]]
  ligand <- c(ligand,temp)
  receptor_tem <- c(receptor_tem,rep(rec,length(temp)))
  cell_type_tem <- c(cell_type_tem,rep(cell,length(temp)))
}
new_non_tad <- data.frame(receptor = receptor_tem,
                              CLIST = rep(0,length(receptor_tem)),
                              MLIST = rep(0,length(receptor_tem)),
                              MLIST.TAD = rep(0,length(receptor_tem)),
                              TAD = rep(0,length(receptor_tem)),
                              VEH = rep(1,length(receptor_tem)),
                              cell_type = cell_type_tem,
                              ligand = ligand)
#saveRDS(new_non_tad,file = 'new_non_tad.rds')










veh <- comp_table[comp_table$VEH==1&comp_table$TAD==0&comp_table$CLIST==0&comp_table$MLIST==0&comp_table$MLISTTAD==0,]
tad <- comp_table[comp_table$VEH==0&comp_table$TAD==1&comp_table$CLIST==0&comp_table$MLIST==0&comp_table$MLISTTAD==0,]
mlist <- comp_table[comp_table$VEH==0&comp_table$TAD==0&comp_table$CLIST==0&comp_table$MLIST==1&comp_table$MLISTTAD==0,]
mlisttad <- comp_table[comp_table$VEH==0&comp_table$TAD==0&comp_table$CLIST==0&comp_table$MLIST==0&comp_table$MLISTTAD==1,]
veh <- veh[,!colnames(veh)%in%c("Freq")]
tad <- tad[,!colnames(tad)%in%c("Freq")]
mlist <- mlist[,!colnames(mlist)%in%c("Freq")]
mlisttad <- mlisttad[,!colnames(mlisttad)%in%c("Freq")]
colnames(veh)[which(names(veh) == "all_int")] <- "receptors"
colnames(tad)[which(names(tad) == "all_int")] <- "receptors"
colnames(mlist)[which(names(mlist) == "all_int")] <- "receptors"
colnames(mlisttad)[which(names(mlisttad) == "all_int")] <- "receptors"

receptor_tem <- c()
cell_type_tem <- c()
ligand <- c()
for (i in 1:nrow(veh)){
  rec <- as.character(veh[i,1])
  cell <- as.character(veh[i,7])
  temp <- dom_veh@linkages$rec_lig[[rec]]
  ligand <- c(ligand,temp)
  receptor_tem <- c(receptor_tem,rep(rec,length(temp)))
  cell_type_tem <- c(cell_type_tem,rep(cell,length(temp)))
}
veh_2 <- data.frame(receptor = receptor_tem,
                       CLIST = rep(0,length(receptor_tem)),
                       MLIST = rep(0,length(receptor_tem)),
                       MLIST.TAD = rep(0,length(receptor_tem)),
                       TAD = rep(0,length(receptor_tem)),
                       VEH = rep(1,length(receptor_tem)),
                       cell_type = cell_type_tem,
                       ligand = ligand)


receptor_tem <- c()
cell_type_tem <- c()
ligand <- c()
for (i in 1:nrow(tad)){
  rec <- as.character(tad[i,1])
  cell <- as.character(tad[i,7])
  temp <- dom_tad@linkages$rec_lig[[rec]]
  ligand <- c(ligand,temp)
  receptor_tem <- c(receptor_tem,rep(rec,length(temp)))
  cell_type_tem <- c(cell_type_tem,rep(cell,length(temp)))
}
tad_2 <- data.frame(receptor = receptor_tem,
                    CLIST = rep(0,length(receptor_tem)),
                    MLIST = rep(0,length(receptor_tem)),
                    MLIST.TAD = rep(0,length(receptor_tem)),
                    TAD = rep(1,length(receptor_tem)),
                    VEH = rep(0,length(receptor_tem)),
                    cell_type = cell_type_tem,
                    ligand = ligand)

receptor_tem <- c()
cell_type_tem <- c()
ligand <- c()
for (i in 1:nrow(mlist)){
  rec <- as.character(mlist[i,1])
  cell <- as.character(mlist[i,7])
  temp <- dom_mlist@linkages$rec_lig[[rec]]
  ligand <- c(ligand,temp)
  receptor_tem <- c(receptor_tem,rep(rec,length(temp)))
  cell_type_tem <- c(cell_type_tem,rep(cell,length(temp)))
}
mlist_2 <- data.frame(receptor = receptor_tem,
                    CLIST = rep(0,length(receptor_tem)),
                    MLIST = rep(1,length(receptor_tem)),
                    MLIST.TAD = rep(0,length(receptor_tem)),
                    TAD = rep(0,length(receptor_tem)),
                    VEH = rep(0,length(receptor_tem)),
                    cell_type = cell_type_tem,
                    ligand = ligand)

receptor_tem <- c()
cell_type_tem <- c()
ligand <- c()
for (i in 1:nrow(mlisttad)){
  rec <- as.character(mlisttad[i,1])
  cell <- as.character(mlisttad[i,7])
  temp <- dom_mlisttad@linkages$rec_lig[[rec]]
  ligand <- c(ligand,temp)
  receptor_tem <- c(receptor_tem,rep(rec,length(temp)))
  cell_type_tem <- c(cell_type_tem,rep(cell,length(temp)))
}
mlisttad_2 <- data.frame(receptor = receptor_tem,
                    CLIST = rep(0,length(receptor_tem)),
                    MLIST = rep(0,length(receptor_tem)),
                    MLIST.TAD = rep(1,length(receptor_tem)),
                    TAD = rep(0,length(receptor_tem)),
                    VEH = rep(0,length(receptor_tem)),
                    cell_type = cell_type_tem,
                    ligand = ligand)

#setwd("E:/WJHLab/domino/local")
#saveRDS(veh_2,file = "veh_2.RDS")
#saveRDS(tad_2,file = "tad_2.RDS")
#saveRDS(mlist_2,file = "mlist_2.RDS")
#saveRDS(mlisttad_2,file = "mlisttad_2.RDS")
#saveRDS(tadmlisttad, file = "tadmlisttad.RDS")
VEH <- readRDS("E:/WJHLab/domino/local/veh_2.RDS")
TAD <- readRDS("E:/WJHLab/domino/local/tad_2.RDS")
MLIST <- readRDS("E:/WJHLab/domino/local/mlist_2.RDS")
MLISTTAD <- readRDS("E:/WJHLab/domino/local/mlisttad_2.RDS")

sobj <- readRDS("E:/WJHLab/singlecell/filtered_sobj_withSampleid.RDS")
unique(sobj@meta.data$sample_id)
sobj$cell_treat <- paste(sobj$cell_type, sobj$group_id, sep = "_")

### All included ligand list

celltype <- c(levels(sobj))
Ligand <- c()
treat <- list(VEH=VEH,TAD=TAD)
cell_treat <- unique(sobj@meta.data$cell_treat)
Ligand1v2 <- data.frame(ligand=NULL,p_val=NULL,avg_log2FC=NULL,pct.1=NULL,pct.2=NULL,p_val_adj=NULL,cell_type=NULL,treatment=NULL,comp=NULL)
for (df_name in names(treat)){
  df <- treat[[df_name]]
  for (i in 1:8){
    temp <- FindMarkers(object = sobj, 
                        ident.1 = paste(celltype[i],"VEH",sep="_"), 
                        ident.2 = paste(celltype[i],"TAD",sep="_"), 
                        features = intersect(rownames(sobj), unique(df$ligand)), 
                        group.by = "cell_treat",
                        min.pct = 0.0001,min.cells.feature = 0, 
                        logfc.threshold = -Inf)
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

celltype <- c(levels(sobj))
Ligand <- c()
treat <- list(TAD=TAD, MLISTTAD=MLISTTAD)
cell_treat <- unique(sobj@meta.data$cell_treat)
LigandTT <- data.frame(ligand=NULL,p_val=NULL,avg_log2FC=NULL,pct.1=NULL,pct.2=NULL,p_val_adj=NULL,cell_type=NULL,treatment=NULL,comp=NULL)
for (df_name in names(treat)){
  df <- treat[[df_name]]
  for (i in 1:8){
    temp <- FindMarkers(object = sobj, ident.1 = paste(celltype[i],"TAD",sep="_"), ident.2 = paste(celltype[i],"MLIST.TAD",sep="_"), features = intersect(rownames(sobj), unique(df$ligand)), group.by = "cell_treat",min.pct = 0.0001,min.cells.feature = 0, logfc.threshold = -Inf)
    temp$cell_type <- rep(celltype[i],nrow(temp))
    temp$treatment <- rep(df_name,nrow(temp))
    Ligand <- c(Ligand,rownames(temp))
    
    LigandTT <- rbind(LigandTT,temp)
  }
}
LigandTT$ligand <- Ligand

new_list1v2 <- sobj@assays$RNA@data[intersect(rownames(Ligand1v2), rownames(sobj@assays$RNA@data)),]
new_list4v5 <- sobj@assays$RNA@data[intersect(rownames(Ligand4v5), rownames(sobj@assays$RNA@data)),]
new_listTT <- sobj@assays$RNA@data[intersect(rownames(LigandTT), rownames(sobj@assays$RNA@data)),]

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

cta <- unique(sobj$cell_treat)[grepl("MLIST.TAD|_[VT]AD$", unique(sobj$cell_treat))]
signalingTT <- matrix(0, ncol = length(cta),
                       nrow = length(rownames(new_listTT)))
rownames(signalingTT) <- rownames(new_listTT)
colnames(signalingTT) <- cta
for(clust in cta){
  sig <- rowMeans(new_listTT[, which(sobj$cell_treat == clust)])
  # ensure the rowmeans order matches the rownames of the signaling matrix
  sig <- sig[rownames(signalingTT)]
  signalingTT[, clust] <- sig
}

setwd("E:/WJHLab/revised_domino/ligandList and signailing")
saveRDS(signaling1v2,"signaling1v2.RDS")
saveRDS(signaling4v5,"signaling4v5.RDS")
saveRDS(signalingTT,"signalingTT.RDS")

signaling1v2 <- readRDS("E:/WJHLab/revised_domino/ligandList and signailing/signaling1v2.RDS")

cl_rec_percent_veh = NULL
ser_receptors = unique(vehTAD_2$receptor)
for(rec in ser_receptors){
  rec_percent <- sapply(
    X = levels(dom_veh@clusters),
    FUN = function(x){
      # percentage of cells in cluster with non-zero expression of receptor
      sum(dom_veh@counts[rec,dom_veh@clusters == x] > 0) / length(dom_veh@counts[rec,dom_veh@clusters == x])
    }
  )
  cl_rec_percent_veh <- rbind(cl_rec_percent_veh, rec_percent)
}
rownames(cl_rec_percent_veh) = ser_receptors

cl_rec_percent_mlist = NULL
ser_receptors = unique(mlistTAD_2$receptor)
for(rec in ser_receptors){
  rec_percent <- sapply(
    X = levels(dom_veh@clusters),
    FUN = function(x){
      # percentage of cells in cluster with non-zero expression of receptor
      sum(dom_mlist@counts[rec,dom_mlist@clusters == x] > 0) / length(dom_mlist@counts[rec,dom_mlist@clusters == x])
    }
  )
  cl_rec_percent_mlist <- rbind(cl_rec_percent_mlist, rec_percent)
}
rownames(cl_rec_percent_mlist) = ser_receptors


library(ggplot2)
heatmap(cl_rec_percent_mlist)





