library(domino)
library(ComplexHeatmap)

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

############################################get 0,1 table results from domino
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

tadlist <- comp_table[,c(1,6,7,8)]
mlistlist <- comp_table[,c(1,4,5,8)]
########################################translating 0,1 into 0,1,2,3 dataframe
vehtad_hm <- data.frame(receptor=NULL,
                        CAF=NULL,
                        DC=NULL,
                        Endothelial=NULL,
                        Granulocyte=NULL,
                        Macrophage=NULL,
                        Myeloid=NULL,
                        Tcell=NULL,
                        Tumor=NULL
)
for (i in 1:length(unique(tadlist$all_int))){
  temp_list <- data.frame(receptor=1,
                          CAF=5,
                          DC=5,
                          Endothelial=5,
                          Granulocyte=5,
                          Macrophage=5,
                          Myeloid=5,
                          Tcell=5,
                          Tumor=5
  )
  rec <- unique(tadlist$all_int)[i]
  temp_list$receptor <- rec
  for (j in 1:length(unique(tadlist$cell_type))){
    ct <- unique(tadlist$cell_type)[j]
    temp <- tadlist[tadlist$all_int==rec&tadlist$cell_type==ct,]
    if (nrow(temp)==0){
      temp_list[[ct]] <- 0
    }  else {
      if (temp$TAD==1&temp$VEH==0) {
        temp_list[[ct]] <- 1
        }  
      if (temp$TAD==0&temp$VEH==1) {
        temp_list[[ct]] <- 2
        }  
      if (temp$TAD==1&temp$VEH==1) {
        temp_list[[ct]] <- 3
        }  
      if (temp$TAD==0&temp$VEH==0) {
        temp_list[[ct]] <- 0
        }
      }
    
  }
  vehtad_hm <- rbind(temp_list,vehtad_hm)
}


mlistmlisttad_hm <- data.frame(receptor=NULL,
                               CAF=NULL,
                               DC=NULL,
                               Endothelial=NULL,
                               Granulocyte=NULL,
                               Macrophage=NULL,
                               Myeloid=NULL,
                               Tcell=NULL,
                               Tumor=NULL
)
for (i in 1:length(unique(mlistlist$all_int))){
  temp_list <- data.frame(receptor=1,
                          CAF=5,
                          DC=5,
                          Endothelial=5,
                          Granulocyte=5,
                          Macrophage=5,
                          Myeloid=5,
                          Tcell=5,
                          Tumor=5
  )
  rec <- unique(mlistlist$all_int)[i]
  temp_list$receptor <- rec
  for (j in 1:length(unique(mlistlist$cell_type))){
    ct <- unique(mlistlist$cell_type)[j]
    temp <- mlistlist[mlistlist$all_int==rec&mlistlist$cell_type==ct,]
    if (nrow(temp)==0){
      temp_list[[ct]] <- 0
    }  else {
      if (temp$MLIST==1&temp$MLISTTAD==0) {
        temp_list[[ct]] <- 1
      }  
      if (temp$MLIST==0&temp$MLISTTAD==1) {
        temp_list[[ct]] <- 2
      }  
      if (temp$MLIST==1&temp$MLISTTAD==1) {
        temp_list[[ct]] <- 3
      }  
      if (temp$MLIST==0&temp$MLISTTAD==0) {
        temp_list[[ct]] <- 0
      }
    }
    
  }
  mlistmlisttad_hm <- rbind(temp_list,mlistmlisttad_hm)
}

write.csv(vehtad_hm,file = "vehtad_hm.csv")
write.csv(mlistmlisttad_hm,file = "mlistmlisttad_hm.csv")





