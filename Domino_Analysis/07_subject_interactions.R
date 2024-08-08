# Differential Interactions
# Calls all interactions identified in the per-patient domino runs and saves them
# in an object class with nested lists of each subject that domino was run on

# DominoInteractions
# Subject
# - Cluster
# - - tfs (active transcription factors)
# - - rec (Receptors with expression significantly correlated with t. factor activity)
# - - lig (Incoming ligands that interact with active receptors)
# - - tfs_rec (linkages of t. factors <- receptors)
# - - rec_lig (linkages of receptors <- ligands)

library(domino)
library(Seurat)
library(dplyr)

sessionInfo()

## Define Domino Interactions Class ###########################################
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
###############################################################################

## Define Functions ###########################################################
# generate the general plots for each domino object
plot_domino <- function(dom, figure_dir, subject = "domino"){
  cell_types <- levels(dom@clusters)
  # incoming signal heatmaps
  pdf(paste0(figure_dir, "/", subject, "_all_incoming_signals_heatmap.pdf"))
  for(cell in cell_types){
    incoming_signaling_heatmap(dom, rec_clust = cell,
                               max_thresh = 2.5)
  }
  dev.off()
  # signaling network
  pdf(paste0(figure_dir, "/", subject, "_signaling_network.pdf"))
  signaling_network(dom, edge_weight = .5, max_thresh = 2.5,
                    main = paste0(subject, " Signaling Network"))
  dev.off()
  # gene netwok
  for(cell in cell_types){
    pdf(paste0(figure_dir, "/", subject, "_domino_", cell, "_gene_network.pdf"),
        onefile=FALSE)
    gene_network(dom, clust = cell, layout = "grid")
    # title(main = paste0(cell, " Gene Network"))
    dev.off()
  }
  # feature heatmap
  pdf(paste0(figure_dir, "/", subject, "_domino_feature_heatmap.pdf"), onefile=FALSE)
  feat_heatmap(dom, norm = TRUE, bool = FALSE)
  dev.off()
  # correlation heatmap
  pdf(paste0(figure_dir, "/", subject, "_domino_correlation_heatmap.pdf"), onefile=FALSE)
  cor_heatmap(dom, bool = FALSE, mark_connections = TRUE)
  dev.off()
}

# write domino interactions object
write_domino_interactions <- 
  function(domino_results, subject_meta, results_names = NULL) {
    # domino results should be a list of domino objects
    # subject_meta is a dataframe that includes the covariates by which the objects could be grouped
    # results_names is what each object should be called, defaults to the first column of subject_meta
    
    if(is.null(results_names)){
      results_names = subject_meta[1,]
    }
    
    subjects = results_names
    subject_features = list() # ligands, receptors, and transcription factors from subject results
    
    for(id in results_names){
      dom <- domino_results[[id]]
      clusters <- levels(dom@clusters)
      
      c_features <- list()
      for(cluster in clusters){
        # empty linkage vectors to fill
        tfs_rec <- c()
        rec_lig <- c()
        
        
        # list of t.factors active in cell type
        tfs <- dom@linkages$clust_tf[[cluster]]
        
        # obtain all receptors linked to active t. factors
        rec <- c()
        for(t in tfs){
          linked_rec <- dom@linkages$tf_rec[[t]]
          for(r in linked_rec){
            # connected t.factor <- receptor are placed next to eachother in vector
            tfs_rec <- c(tfs_rec, t, r)
          }
          rec <- c(rec, linked_rec)
        }
        
        # limit to unique entries
        tfs <- unique(tfs)
        rec <- unique(rec)
        
        # obtain all incoming ligands that interact with the receptors
        # limited to those present in data set
        allowed_lig <- rownames(dom@cl_signaling_matrices[[cluster]])
        
        lig <- c()
        for(r in rec){
          linked_lig <- dom@linkages$rec_lig[[r]]
          for(l in linked_lig){
            if(length(which(allowed_lig == l))){
              rec_lig = c(rec_lig, r, l)
              lig = c(lig, l)
              }
            }
          }
        # limit to unique ligands
        lig <- unique(lig)
        
        # remove the elipses from the t.factor names
        tfs <- gsub("\\.\\.\\.", "", tfs)
        tfs_rec <- gsub("\\.\\.\\.", "", tfs_rec)
        
        # stitch linked t.factors-receptors, receptors-ligands
        int_tfs_rec <- c()
        int_rec_lig <- c()
        for(i in 0:((length(tfs_rec)/2)-1)){
          # count by twos and paste together with a <- denoting direction
          s <- i * 2
          interact <- paste(tfs_rec[1 + s], tfs_rec[2 + s], sep = " <- ")
          int_tfs_rec <- c(int_tfs_rec, interact)
        }
        for(i in 0:((length(rec_lig)/2)-1)){
          # count by twos and paste together with a "<-" denoting direction
          s <- i * 2
          interact <- paste(rec_lig[1 + s], rec_lig[2 + s], sep = " <- ")
          int_rec_lig <- c(int_rec_lig, interact)
        }
        
        # save the features of this cluster
        c_features[[cluster]] <- list(
          "tfs" = tfs,
          "rec" = rec,
          "lig" = lig,
          "tfs_rec" = int_tfs_rec,
          "rec_lig" = int_rec_lig
        )
        
        }
      subject_features[[id]] <- c_features
    }
    subjects <- factor(subjects)
    
    return(
      DominoInteractions(
        subjects = subjects,
        subject_meta = subject_meta,
        subject_features = subject_features
      )
    )
}

# count occurrences of interactions among subjects in each group
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

###############################################################################

# results directory
result_dir <- "processed_data/07_tadalafil_interactions"
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}
writeLines(capture.output(sessionInfo()),
           con = paste0(result_dir, "/session_info.txt"))

# figure directory
figure_dir <- "figures/07_tadalafil_interactions"
if(!dir.exists(figure_dir)){
  dir.create(figure_dir)
}

# data frame of patient info annotations
pt_meta <- data.frame(
  "Treatment" = c("CLIST","MLIST","MLIST.TAD","TAD","VEH")
)

write.csv(pt_meta, file = paste0(result_dir, "/tadalafil_meta.csv"))

# load domino results for patients of interest
domino_res_dir <- "processed_data/06_domino"
domino_res_files <- list.files(domino_res_dir)

domino_res <- list()
for(i in domino_res_files){
  p <- gsub("", "",
            gsub("_domino_unbuilt\\.rds$", "", i))
  domino_res[[p]] <-
    readRDS(
      paste0(domino_res_dir, "/", i)
      )
  # domino object build parameters
  domino_res[[p]] <-
    build_domino(
      dom = domino_res[[p]],
      min_tf_pval = .001,
      max_tf_per_clust = Inf,
      max_rec_per_tf = Inf,
      rec_tf_cor_threshold = .25
    )
}

a <- 1
for(dom in domino_res){
  dom_name <- names(domino_res)[a]
  cat(dom_name)
  plot_domino(dom, 
              figure_dir = figure_dir,
              subject = dom_name)
  
  a <- a + 1
}
a <- 1

dom_diff <- write_domino_interactions(
  domino_results = domino_res,
  subject_meta = tadalafil_meta,
  results_names = tadalafil_meta$Treatment
)

# save result of interactions list
saveRDS(dom_diff,
        file = paste0(result_dir, "/tadalafil_domino_features.rds"))

# clear global environment to prevent the workspace from being saved
rm(list = ls())

