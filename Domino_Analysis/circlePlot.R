library(domino)

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
setwd("E:/WJHLab/revised_domino/circusPlot/circusPlot2")
domino_MLIST <- readRDS(file = "E:/WJHLab/revised_domino/new_unbuilt/MLIST_domino_unbuilt.rds")
domino_MLISTTAD <- readRDS(file = "E:/WJHLab/revised_domino/new_unbuilt/MLISTTAD_domino_unbuilt.rds")
domino_TAD <- readRDS(file = "E:/WJHLab/revised_domino/new_unbuilt/TAD_domino_unbuilt.rds")
domino_VEH <- readRDS(file = "E:/WJHLab/revised_domino/new_unbuilt/VEH_domino_unbuilt.rds")

domino_TAD <-
  build_domino(
    dom = domino_TAD,
    min_tf_pval = .001,
    max_tf_per_clust = Inf,
    max_rec_per_tf = Inf,
    rec_tf_cor_threshold = .25
  )
domino_VEH <-
  build_domino(
    dom = domino_VEH,
    min_tf_pval = .001,
    max_tf_per_clust = Inf,
    max_rec_per_tf = Inf,
    rec_tf_cor_threshold = .25
  )
domino_MLIST <-
  build_domino(
    dom = domino_MLIST,
    min_tf_pval = .001,
    max_tf_per_clust = Inf,
    max_rec_per_tf = Inf,
    rec_tf_cor_threshold = .25
  )
domino_MLISTTAD <-
  build_domino(
    dom = domino_MLISTTAD,
    min_tf_pval = .001,
    max_tf_per_clust = Inf,
    max_rec_per_tf = Inf,
    rec_tf_cor_threshold = .25
  )

active_rec <- lapply(
  domino_TAD@linkages$clust_tf,
  function(x){
    linked = unlist(domino_TAD@linkages$tf_rec[names(domino_TAD@linkages$tf_rec) %in% x])
    return(unique(linked))
  }
)

ggplot_col_gen = function(n){
  hues = seq(15, 375, length = n + 1)
  return(hcl(h = hues, l = 65, c = 100)[1:n])
}

circos_ligand_receptor = function(dom, receptor, ligand_expression_threshold = 0.01, cell_idents = NULL, cell_colors = NULL){
  require(circlize)
  require(ComplexHeatmap)
  
  ligands <- dom@linkages$rec_lig[[receptor]]
  ligands <- ligands[ligands %in% rownames(dom@z_scores)]
  if(length(ligands) == 0){stop(paste0("no ligands of ", receptor, " are present in the expression data"))}
  signaling_df <- NULL
  
  if(is.null(cell_idents)){
    # default to all cluster labels in domino object in alphabetical order
    cell_idents <- sort(unique(dom@clusters))
  }
  
  # obtain expression values from cl_signaling matrices
  active_chk <- sapply(
    dom@linkages$clust_tf,
    function(x){
      linked = unlist(dom@linkages$tf_rec[names(dom@linkages$tf_rec) %in% x])
      return(receptor %in% linked)
    }
  )
  if(sum(active_chk)){
    # obtain a signaling matrix where receptor is active
    active_cell <- names(active_chk[active_chk == TRUE])
    sig <- dom@cl_signaling_matrices[active_cell][[1]]
    cell_names <- gsub("^L_", "", colnames(sig))
    for(l in ligands){
      df <- data.frame(
        "origin" = paste0(cell_names, "-", l),
        "destination" = receptor,
        "mean.expression" = unname(sig[rownames(sig) == l,])
      )
      signaling_df <- rbind(signaling_df, df)
    }
  } else {stop(paste0("No clusters have active ", receptor, " signaling"))}
  
  signaling_df$mean.expression[is.na(signaling_df$mean.expression)] <- 0
  # create a scaled mean expression plot for coord widths greater than 1
  # by dividing by the max expression [range (0-1)]
  # scaled.mean will only be used when the max expression is > 1
  signaling_df$scaled.mean.expression <- signaling_df$mean.expression/max(signaling_df$mean.expression)
  
  # exit function if no ligands are expressed above ligand expression threshold
  if(sum(signaling_df[["mean.expression"]] > ligand_expression_threshold) == 0){
    stop(paste0("No ligands of ", receptor, " exceed ligand expression threshold."))
  }
  
  # initialize chord diagram with even ligand arcs
  arc_df <- signaling_df[, c("origin", "destination")]
  arc_df["ligand.arc"] <- 1
  # receptor arc will always sum to 4 no matter how many ligands and cell idents are plotted
  arc_df["receptor.arc"] <- 4 / (nrow(signaling_df)) 
  
  # name grouping based on [cell_ident]
  nm <- c(receptor, arc_df$origin)
  group <- structure(c(nm[1], gsub("-.*", "", nm[-1])), names = nm)
  
  # order group as a factor with the receptor coming first
  group <- factor(group,
                  levels = c(receptor,
                             sort(unique(gsub("-.*", "", nm))[-1]) # alphabetical order of the other cell idents
                  ))
  
  # colors for ligand chords
  lig_colors <- ggplot_col_gen(length(ligands))
  names(lig_colors) <- ligands
  
  # colors for [cell_ident] arcs
  if(is.null(cell_colors)){
    cell_colors <- ggplot_col_gen(length(cell_idents))
    names(cell_colors) <- cell_idents
  }
  
  grid_col <-  c("#FFFFFF") # hide the arc corresponding to the receptor by coloring white 
  for(i in 1:length(ligands)){
    grid_col <- c(grid_col, rep(lig_colors[i], length(cell_idents)))
  }
  names(grid_col) <- c(receptor, signaling_df$origin)
  circos.clear()
  circos.par(start.degree = 0)
  
  chordDiagram(arc_df, group = group, grid.col = grid_col,
               link.visible = FALSE, # hide default chords
               annotationTrack = c("grid"),
               preAllocateTracks = list(
                 track.height = mm_h(4),
                 track.margin = c(mm_h(2), 0)),
               big.gap = 2
  )
  
  for(send in signaling_df$origin){
    if(signaling_df[signaling_df$origin == send,][["mean.expression"]] > ligand_expression_threshold){
      if(max(signaling_df[["mean.expression"]]) > 1){
        expr <- signaling_df[signaling_df$origin == send,][["scaled.mean.expression"]]
        max_width <- signif(max(signaling_df[["mean.expression"]]), 2)
      } else {
        expr <- signaling_df[signaling_df$origin == send,][["mean.expression"]]
        max_width <- 1
      }
      
      circos.link(send, 
                  c(0.5 - (expr/2), 0.5 + (expr/2)), 
                  receptor, 2,
                  col = paste0(grid_col[[send]], "88"))
    }
  }
  sector_names <- get.all.sector.index()
  cell_sectors <- cell_idents[cell_idents %in% gsub("-.*", "", sector_names)]
  
  for(cell in cell_sectors){
    row_pick <- sector_names[grepl(paste0("^", cell), sector_names)]
    
    if(length(row_pick)){
      highlight.sector(
        sector_names[grepl(paste0("^", cell, "-"), sector_names)],
        track.index = 1, col = cell_colors[[cell]],
        text = cell,
        cex = 1, facing = "inside",
        text.col = "black", niceFacing = FALSE,
        text.vjust = -1.5
      )
    }
  }
  # highlight receptor sector
  highlight.sector(
    sector_names[grepl(paste0("^", receptor, "$"), sector_names)], 
    track.index = 1, col = "#FFFFFF",
    text = receptor, cex = 1.5, facing = "clockwise",
    text.col = "black", niceFacing = TRUE,
    pos = 4
  )
  # create legends
  lgd_cells = Legend(
    at = as.character(cell_idents), type = "grid", 
    legend_gp = gpar(fill = cell_colors),
    title_position = "topleft", title = "cell identity"
  )
  lgd_ligands = Legend(
    at = ligands, type = "grid", 
    legend_gp = gpar(fill = lig_colors),
    title_position = "topleft", title = "ligand"
  )
  chord_width <- 10/(4 + length(cell_idents)*length(ligands))
  lgd_chord = Legend(
    at = c(ligand_expression_threshold, max_width), 
    col_fun = colorRamp2(c(ligand_expression_threshold, max_width), c("#DDDDDD", "#DDDDDD")),
    legend_height = unit(chord_width, "in"),
    title_position = "topleft", title = "ligand expression"
  )
  lgd_list_vertical = packLegend(lgd_cells, lgd_ligands, lgd_chord)
  draw(lgd_list_vertical, 
        x = unit(0.02, "npc"), y = unit(0.98, "npc"),
       just = c("left", "top"))
}

vehlist <- c("Csf1r","Hcst","Mrc1","Notch1","Notch2","Sirpa","Pirb","Tgfbr1","Adora2a","Ccr2","Ccr5","Il10rb","Il17ra","I4ra","Il6ra","Tgfbr2")
tadlist <- c("Plxnb2","Clec2d","Ltbr","Tnfrsf23","Tnfrsf4","Lifr")
mlistlist <- c("Notch1","Notch2")
mlisttadlist <- c("Tnfrsf23","Cd27","Cd96","Cxcr3","Il12rb2","Pvr","Tnfrsf4","Lilra6")
length(mlisttadlist)
for (receptor in mlisttadlist){
  pdf(paste0(receptor, "_MLISTTAD.pdf"),height = 22, width = 22)
  circos_ligand_receptor(domino_MLISTTAD,receptor,ligand_expression_threshold = 0,
                         cell_idents = NULL, cell_colors = NULL)
  dev.off()
}

sum(sapply(domino_TAD@linkages$tf_rec, function(x) "Insr" %in% x))
pdf("Fcgrt_TAD.pdf",height = 22, width = 22)
circos_ligand_receptor(domino_TAD,"Fcgrt",ligand_expression_threshold = 0,
                       cell_idents = NULL, cell_colors = NULL)
dev.off()


