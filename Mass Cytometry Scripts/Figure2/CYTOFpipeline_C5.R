#Nicole Gross
#R version 4.0.2 (2020-06-22)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 18.04.5 LTS

rm(list = ls())


####READ and CLUSTER FUNCTIONS####
returnfcs <- function(FDR_cutoff=.05,
                      metaDataFile='~/',
                      panelDataFile='~/',
                      dataDirectory='~/',
                      shape_conditions=NULL,
                      color_conditions=NULL){
  #This function generates an fcs file, subtype_markers, colors and shapes for clustering 
  require(scales);require(readxl);require(dplyr);require(flowCore)
  ##directory and metadatafile checking
  if(!dir.exists(dataDirectory)) {stop('ERR: cannot find data directory')}
  if(!file.exists(metaDataFile)) {stop('ERR: cannot find metadata.xlsx or .csv file')}
  ##readin metadata and clean
  ifelse(grepl(metaDataFile,pattern='.xls'),md <- read_excel(metaDataFile),md <- read.csv(metaDataFile,header = TRUE))#must be in xl format or csv
  md$condition <- factor(md$condition)
  md$batch <- factor(md$batch)
  md$run <- factor(md$run)
  rownames(md) = md$sample_id;md$sample_id <- md$sample_id
  #Make sure all files in metadata present in datadirectory
  if(!all(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])){
    print(paste('ERR: not all filenames in metadata present in data folder - missing',md$file_name[!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])],'Subsetting...'))
    md <- md[-c(!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])),]
  }
  ##Define shapes for conditions
  if(is.null(shape_conditions)){shape_conditions <- c(0:25)[1:length(levels(md$condition))]}#can specify as long as number is same
  if(length(shape_conditions)!=length(levels(md$condition))){stop(paste0('ERR no. shapes specified is less than no. of conditions (',length(levels(md$condition)),')'))}
  names(shape_conditions) <- levels(md$condition)
  ## Define colors for the conditions
  if(is.null(color_conditions)){color_conditions <- hue_pal()(length(levels(md$condition)))}#can specify as long as number is same
  if(length(color_conditions)!=length(levels(md$condition))){stop(paste0('ERR no. shapes specified is less than no. of conditions (',length(levels(md$condition)),')'))}
  ## read fcs
  fcs_raw <- read.flowSet(md$file_name, path = dataDirectory, transformation = FALSE, truncate_max_range = FALSE)
  
  panel <- read_excel(panelDataFile)
  head(data.frame(panel))
  ## Replace problematic characters
  panel$Metal <- gsub('-', '_', panel$Metal)
  panel_fcs <- pData(parameters(fcs_raw[[1]]))
  panel_fcs$desc <- gsub('-', '_', panel_fcs$desc)
  panel_fcs$desc[is.na(panel_fcs$desc)] <- paste0('NA_',which(is.na(panel_fcs$desc))) #was labelled 'todo'(mapping based on isotope for now just getting rid of NA keeping rownum) 
  # use panel$Antigen to fix description in panel_fcs
  # use metal+isotope as mapping between panel from xlsx and panel from the fcs files
  rownames(panel_fcs) = panel_fcs$name
  panel_fcs[paste0(panel$Metal,panel$Isotope,'Di'),2] <- panel$Antigen
  ## Replace paramater data in flowSet
  pData(parameters(fcs_raw[[1]])) <- panel_fcs
  ## Define variables indicating marker types
  subtype_markers <- panel$Antigen[panel$Subtype == 1]
  functional_markers <- panel$Antigen[panel$Functional == 1]
  if(!all(subtype_markers %in% panel_fcs$desc)){stop('ERR: Not all subtype_markers in panel_fcs$desc (isotopes)')}
  if(!all(functional_markers %in% panel_fcs$desc)){stop('ERR: Not all functional_markers in panel_fcs$desc (isotopes)')}
  ## arcsinh transformation and column subsetting
  fcs <- fsApply(fcs_raw, function(x, cofactor = 5){
    colnames(x) <- panel_fcs$desc
    expr <- exprs(x)
    expr <- asinh(expr[, union(subtype_markers,functional_markers)] / cofactor)
    expr[!is.finite(expr)] <- NA #convert inf to NA
    expr<-na.omit(expr) #remove NA
    exprs(x) <- expr
    x
  })
  
  sample_ids <- rep(md$sample_id, fsApply(fcs, nrow))
  return(list('fcs'=fcs,
              'subtype_markers'=subtype_markers,
              'functional_markers'=functional_markers,
              'shape_conditions'=shape_conditions,
              'color_conditions'=color_conditions,
              'sample_ids'=sample_ids,
              'meta_data'=md))
}

clusterfcs <- function(fcs=output$fcs,
                       subtype_markers = output$subtype_markers,
                       seed=1234,plottitle='consensus_plots',
                       numclusters=40){
  ## Cell population identification with FlowSOM and ConsensusClusterPlus
  require(dplyr);require(FlowSOM);require(ConsensusClusterPlus)
  set.seed(seed)
  som <- ReadInput(fcs,transform = FALSE, scale = FALSE) %>% BuildSOM(colsToUse = subtype_markers)
  ## Get the cell clustering into 100 SOM codes
  cell_clustering_som <- som$map$mapping[,1]
  ## Metaclustering into numclusters with ConsensusClusterPlus
  codes <- som$map$codes
  mc <- ConsensusClusterPlus(t(codes), maxK = numclusters, reps = 100,
                             pItem = 0.9, pFeature = 1, title = plottitle, 
                             plot = "png", clusterAlg = "hc", 
                             innerLinkage = "average", finalLinkage = "average",
                             distance = "euclidean", seed = 1234)
  
  ## Get cluster ids for each cell
  code_clustering <- mc[[numclusters]]$consensusClass#metaclusters consensus
  cell_clustering <- code_clustering[cell_clustering_som]#cell clustering from som
  return(list('code_clustering'=code_clustering,'cell_clustering'=cell_clustering,'metaclusters'=mc))
}


####CLUSTER HEATMAP FUNCTIONS ####
plot_clustering_heatmap_wrapper <- function(fcs, cell_clustering, nclusters=40,
                                            color_clusters=clustercolors, cluster_merging = NULL, 
                                            subtype_markers,
                                            clusterMergeFile=NULL,
                                            fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  #if((color_clusters)=='auto'){color_clusters <- hue_pal()(nclusters)}
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  ## Calculate the mean expression##################################################
  pdf(fileName, width=8, height=11) 
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  color_heat <- colorRampPalette(brewer.pal(n = 9, name = "YlOrBr"))(100)
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  labels_row <- paste0(expr01_mean$cell_clustering, " (", clustering_prop ,
                       "%)")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  ## Annotation for the merged clusters
  
  if(!is.null(clusterMergeFile)){
    ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Merged <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Merged <- color_clusters2
  }
  
  p <- pheatmap(expr_heat, color = rev(colorRampPalette(brewer.pal(n = 11, name = "RdYlBu"))(100)), 
                cluster_cols = T,
                cluster_rows = T, 
                labels_row = labels_row,
                #scale="column",
                display_numbers = FALSE, number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 8,
                cellheight = 8
  )
  dev.off() 
  print('Colors:')
  print(color_clusters)
  print(p);return(p)
}

plot_clustering_heatmap_wrapper2 <- function(fcs, cell_clustering, nclusters=40,
                                             color_clusters=clustercolors,
                                             subtype_markers,
                                             fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  #if((color_clusters)=='auto'){color_clusters <- hue_pal()(nclusters)}
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  ## Calculate the mean expression##################################################
  pdf(fileName, width=8, height=11) 
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  #labels_row <- expr01_mean$cell_clustering
  
  labels_row <- paste0(expr01_mean$cell_clustering, " ")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  p <- pheatmap(expr_heat, 
                color = c(rep(magma(100)[1],25),magma(100)[1:100]), 
                cluster_cols = T,
                cluster_rows = F, 
                labels_row = labels_row,
                #scale="column",
                display_numbers = F, 
                number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 8,
                cellheight = 8,
                border_color = "black",
                annotation_legend = F
  )
  dev.off() 
  print('Colors:')
  print(color_clusters)
  print(p);return(p)
}

####DIAGNOSTICS####
makeDiagnosticPlots = function(exprData, 
                               md = output$meta_data,
                               sample_ids = output$sample_ids,
                               fcs = output$fcs,
                               subtype_markers = output$subtype_markers,
                               color_conditions = clustercolors,
                               shape_conditions = c(1:13),
                               fileName = 'diagnostics.pdf', 
                               tit = '', 
                               fun = mean)
{
  pdf(file = fileName)
  
  # plot 1
  ggdf <- data.frame(sample_id = sample_ids, exprData)
  ggdf <- melt(ggdf, id.var = 'sample_id', value.name = 'expression', 
               variable.name = 'antigen')
  mm <- match(ggdf$sample_id, md$sample_id)
  ggdf$condition <- md$condition[mm]
  print(ggplot(ggdf, aes(x = expression, color = condition, group = sample_id)) + 
          geom_density() +
          facet_wrap(~ antigen, nrow = 4, scales = 'free') + theme_bw() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                strip.text = element_text(size = 7),
                axis.text = element_text(size = 5)) + 
          scale_color_manual(values = color_conditions) )
  # plot 2
  
  ## Spot check - number of cells per sample
  cell_table <- table(sample_ids)
  ggdf <- data.frame(sample_id = names(cell_table), 
                     cell_counts = as.numeric(cell_table))
  mm <- match(ggdf$sample_id, md$sample_id)
  ggdf$condition <- md$condition[mm]
  print(ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = condition)) + 
          geom_bar(stat = 'identity') + 
          geom_text(aes(label = cell_counts), hjust = 0.5, vjust = -0.5, size = 2.5) + 
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  
          scale_fill_manual(values = color_conditions, drop = FALSE) + 
          scale_x_discrete(drop = FALSE))
  
  ## Multi-dimensional scaling plot to show similarities between samples
  # plot 3
  
  ## Get the mean marker expression per sample##############################################
  expr_mean_sample_tbl <- data.frame(sample_id = sample_ids, exprData) %>%
    group_by(sample_id) %>%  summarize_all(funs(fun))
  expr_mean_sample <- t(expr_mean_sample_tbl[, -1])
  colnames(expr_mean_sample) <- expr_mean_sample_tbl$sample_id
  mds <- plotMDS(expr_mean_sample, plot = FALSE)
  ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
                     sample_id = colnames(expr_mean_sample))
  mm <- match(ggdf$sample_id, md$sample_id)
  ggdf$condition <- md$condition[mm]
  print(ggplot(ggdf, aes(x = MDS1, y = MDS2, color = condition, shape = condition)) +
          geom_point(size = 2, alpha = 0.8) +
          geom_label_repel(aes(label = sample_id)) +
          theme_bw() +
          scale_color_manual(values = color_conditions) +
          scale_shape_manual(values = shape_conditions) +
          coord_fixed())
  
  # plot 4
  ## Can see differences between tissues as well as conditions
  ## Column annotation for the heatmap
  mm <- match(colnames(expr_mean_sample), md$sample_id)
  annotation_col <- data.frame(condition = output$meta_data$condition[mm], 
                               row.names = colnames(expr_mean_sample))
  annotation_colors <- list(condition = color_conditions[1:length(levels(annotation_col$condition))])
  
  ## Colors for the heatmap
  color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
  pheatmap(expr_mean_sample, color = color, display_numbers = TRUE,
           number_color = "black", fontsize_number = 5, 
           annotation_col = levels(annotation_col), main = tit,
           annotation_colors = annotation_colors, clustering_method = "average")
  
  # plot 5
  
  ## Define a function that calculates the Non-Redundancy Score per sample
  NRS <- function(x, ncomp = 3){
    pr <- prcomp(x, center = TRUE, scale. = FALSE)
    score <- rowSums(outer(rep(1, ncol(x)), pr$sdev[1:ncomp]^2) * 
                       abs(pr$rotation[,1:ncomp]))
    return(score)
  }
  
  ## Calculate the score
  ## May want to do the same with other markers
  nrs_sample <- fsApply(fcs[, subtype_markers], NRS, use.exprs = TRUE)
  rownames(nrs_sample) <- md$sample_id
  nrs <- colMeans(nrs_sample, na.rm = TRUE)
  
  ## Plot the NRS for ordered markers
  ## May be helpful to look at tissue instead of condition
  subtype_markers_ord <- names(sort(nrs, decreasing = TRUE))
  nrs_sample <- data.frame(nrs_sample)
  nrs_sample$sample_id <- rownames(nrs_sample)
  ggdf <- melt(nrs_sample, id.var = "sample_id",
               value.name = "nrs", variable.name = "antigen")
  ggdf$antigen <- factor(ggdf$antigen, levels = subtype_markers_ord)
  mm <- match(ggdf$sample_id, md$sample_id)
  ggdf$condition <- md$condition[mm]
  print(ggplot(ggdf, aes(x = antigen, y = nrs)) +
          geom_point(aes(color = condition), alpha = 0.9,
                     position = position_jitter(width = 0.3, height = 0)) +
          geom_boxplot(outlier.color = NA, fill = NA) +
          stat_summary(fun.y = "mean", geom = "point", shape = 21, fill = "white") +
          theme_bw() + ggtitle(tit)+ 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) )
  # scale_color_manual(values = color_conditions)
  
  dev.off()
}


####CLUSTER HISTO####

plot_clustering_distr_wrapper <- function(expr = expr, 
                                          cell_clustering){
  # Calculate the median expression
  cell_clustering <- factor(cell_clustering)
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  # Calculate cluster frequencies
  freq_clust <- table(cell_clustering)
  freq_clust <- round(as.numeric(freq_clust)/sum(freq_clust)*100, 2)
  cell_clustering <- factor(cell_clustering,
                            labels = levels(cell_clustering))
  ### Data organized per cluster
  ggd <- melt(data.frame(cluster = cell_clustering, expr),
              id.vars = "cluster", value.name = "expression",
              variable.name = "antigen")
  ggd$antigen <- factor(ggd$antigen, levels = colnames(expr))
  ggd$reference <- "no"
  ### The reference data
  ggd_bg <- ggd
  ggd_bg$cluster <- "reference"
  ggd_bg$reference <- "yes"
  
  ggd_plot <- rbind(ggd, ggd_bg)
  ggd_plot$cluster <- factor(ggd_plot$cluster,
                             levels = c(levels(cell_clustering)[rev(cluster_rows$order)], "reference"))
  
  ggplot() +
    geom_density_ridges(data = ggd_plot, aes(x = expression, y = cluster,
                                             color = reference, fill = reference), alpha = 0.3) +
    facet_wrap( ~ antigen, scales = "free_x", nrow = 2) +
    theme_ridges() +
    theme(axis.text = element_text(size = 7),  
          strip.text = element_text(size = 7), legend.position = "none")
  
}
####UMAP####
do_umap <- function(fcs,subtype_markers,sample_ids,cell_clustering,metadata,
                    clusterMergeFile='~/Desktop/Analysis BETi Run2/Config/BETi_merged.xlsx',
                    seed = 1234, ncells=2000,sample_subset=NULL){
  require(umap);require(flowCore);require(readxl)
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Create vector to later find and skip duplicates
  dups <- duplicated(expr[, subtype_markers])
  dups <- which(!(dups))## Find and skip duplicates
  ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
  ## New clustering1m
  mm <- match(cell_clustering, cluster_merging$original_cluster)
  cell_clustering1m <- cluster_merging$new_cluster[mm]
  ## Create a data frame of sample_ids and cell_clustering1m
  dtf<-data.frame(ids=sample_ids,type=cell_clustering1m)
 
  ## Data subsampling: create indices by sample
  inds <- split(1:length(sample_ids), sample_ids) #to get original indexes belonging to each cluster
  samplenames <- names(inds) #create a name vector of the files

  custom.settings = umap.defaults
  custom.settings$seed = seed
 
  ifelse(is.null(sample_subset),
         umap_ncells <- pmin(table(sample_ids), ncells),
         umap_ncells <- pmin(table(sample_ids), ncells)[sample_subset]
  )
  if(!is.null(sample_subset)){inds <- inds[sample_subset]}
  umap_inds <- lapply(names(inds), function(i){
    s <- sample(inds[[i]], umap_ncells[i], replace = FALSE)
    intersect(s, dups)
  })
  set.seed(seed)
  umap_inds <- unlist(umap_inds)
  umap_out <- umap(expr[umap_inds, subtype_markers], config = custom.settings, method = 'naive')
  umapRes2D = data.frame(umap1 = umap_out$layout[, 1], umap2 = umap_out$layout[, 2], 
                         expr[umap_inds, subtype_markers],
                         sample_id = sample_ids[umap_inds], cell_clustering = factor(cell_clustering1m[umap_inds]), check.names = FALSE)
  
  return(umapRes2D)
}

plotUmap <- function(umapRes,seed=1234,neighbors=10,midpoint,color_clusters=clustercolors,code_clustering,subtype_markers=NULL)
{require(umap);require(ggplot2);require(viridis);require(ggrepel)
  #if((color_clusters)=='auto'){color_clusters <- hue_pal()(length(unique(code_clustering)))}
  custom.settings = umap.defaults
  custom.settings$seed = seed
  custom.settings$n.neighbors = neighbors
  ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = cell_clustering)) +
    geom_point(size = 1) +
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    
    scale_color_manual(values = color_clusters, name="CLUSTERS") +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  
  print(ggp)
  #other options
  print(ggp + facet_wrap(~ condition)+ggtitle('CONDITIONS'))
  print(ggp + facet_wrap(~ batch)+ggtitle('BATCH'))
  print(ggp + facet_wrap(~ sample_id, ncol = 3)+ggtitle('SAMPLE'))
  ggp2 <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = sample_id)) +
    geom_point(size = 1) +
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    
    scale_color_manual(values = color_clusters, name="CLUSTERS") +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  
  print(ggp2)
  #can specify which markers to display
  if(!is.null(subtype_markers)){
    for(i in subtype_markers)
    {
      ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = umapRes[,i])) +
        geom_point(size = 1) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        scale_color_gradient2(i, low="dark blue",mid="white",high="dark red", midpoint = mean(unlist(umapRes[,i])))
      print(ggp)
    }
  }
}

####OTHER PACKAGES####
library(reshape2)
library(randomcoloR)
library(pals)
library(ggplot2)
library(Hmisc) 
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(corrplot)
library(readxl)
library(ggridges)
library(dplyr)
library(scales)
library(flowCore)
library(FlowSOM)
library(ConsensusClusterPlus)
library(umap)
library(limma)
library(DelayedMatrixStats)
#======================
#     RUNNING DATA
#======================


####DATA LOADING AND CLUSTERING####

#set working directory based on path of script file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
workd<-getwd()

###Start here if loading output rds

#load up saved output RDS
output<-readRDS('backup_output.rds')

###Skip lines 543-553 if loading output rds
#read and cluster
output <- returnfcs(metaDataFile = paste0(workd,"/config/","metadata.xlsx"),
                    panelDataFile = paste0(workd,"/config/","panel.xlsx"),
                    dataDirectory = paste0(workd,'/data'))

output[8:10] <- clusterfcs(numclusters=35)

names(output)[8:10] <- c('code_clustering','cell_clustering','metaclusters')

#assign clusters
clusterMergeFile = paste0(workd,"/config/",'merged.xlsx')
cluster_merging <- read_xlsx(clusterMergeFile)

#set up factor levels

samplevels<-c("C_1",
              "C_2",
              "I_2",
              "I_3",
              "IT_3",
              "IT_4",
              "L_4",
              "L_5",
              "LI_5",
              "LI_6",
              "C_3",
              "C_4",
              "I_4",
              "I_5",
              "IT_5",
              "IT_6",
              "L_6",
              "L_7",
              "LI_7",
              "LIT_1",
              "C_5",
              "C_6",
              "I_6",
              "I_7",
              "IT_7",
              "L_1",
              "LI_1",
              "LI_2",
              "LIT_2",
              "LIT_3",
              "C_7",
              "I_1",
              "IT_1",
              "IT_2",
              "L_2",
              "L_3",
              "LI_3",
              "LI_4",
              "LIT_4",
              "LIT_5",
              "C_8",
              "C_13",
              "I_11",
              "IT_9",
              "IT_14",
              "L_12",
              "LI_10",
              "LI_14",
              "LIT_11",
              "C_9",
              "C_14",
              "I_12",
              "IT_10",
              "L_8",
              "L_13",
              "LI_11",
              "IT_13",
              "LIT_12",
              "C_10",
              "I_8",
              "I_13",
              "IT_11",
              "L_9",
              "L_14",
              "LI_12",
              "LIT_9",
              "LIT_13",
              "C_11",
              "I_9",
              "I_14",
              "IT_12",
              "L_10",
              "LI_8",
              "LI_13",
              "LIT_10",
              "LIT_14",
              "C_12",
              "I_10",
              "IT_8",
              "LIT_8",
              "L_11",
              "LI_9")

grouplevels=c("C",
              "I",
              "IT",
              "L",
              "LI",
              "LIT")

#clusterlevels=c(1:35) #assigning cluster levels before annotation
clusterlevels = c("B",
                  "Th_I",
                  "Th_II",
                  "Tc",
                  "DC_I",
                  "DC_II",
                  "NK",
                  "TAM1",
                  "TAM2_I",
                  "TAM2_II",
                  "Gran_I",
                  "Gran_II",
                  "MMDSC",
                  "Endothelial",
                  "Tumor",
                  "UA")

runlevels= c("1","2")
                  
clustercolors <- as.character(c(cols25(n=25),alphabet(n=26),alphabet2(n=26),glasbey(n=32)))

###Skip lines 670-672 if loading output rds
mm1 <- match(output$cell_clustering, cluster_merging$original_cluster)
cell_clustering1m <- cluster_merging$new_cluster[mm1]
output$cell_clustering1m <- cell_clustering1m

#metacluster heatmap
plot_clustering_heatmap_wrapper2(fcs=output$fcs,
                                 color_clusters = clustercolors,
                                 cell_clustering = factor(output$cell_clustering1m, levels=clusterlevels), 
                                 subtype_markers=output$subtype_markers,
                                 fileName = 'clusteringheatmap_final.pdf');dev.off()


#save output
#saveRDS(output, file = "backup_output.rds")


expr <- fsApply(output$fcs, exprs) #create expression matrix
rng <- colQuantiles(expr, probs = c(0.01, 0.99))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1])) #scaling 0-1
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1 # need this expr object later 

## Umap

#to load umap rds
umapRes<-readRDS('backup_umap.rds')

###skip lines 698-708 if loading UMAP rds
umapRes <- do_umap(fcs=output$fcs,subtype_markers = output$subtype_markers,
                   sample_ids = output$sample_ids,cell_clustering = output$cell_clustering, metadata=output$metadata,
                   clusterMergeFile=clusterMergeFile,
                   seed = 1234, ncells=500,sample_subset=NULL)

mm <- match(as.character(umapRes$sample_id), as.character(output[["meta_data"]]$sample_id))
umapRes$condition <- factor(output[["meta_data"]]$condition[mm], levels=grouplevels)
umapRes$batch <- output[["meta_data"]]$batch[mm]
umapRes$run <- output[["meta_data"]]$run[mm]
umapRes$sample_id <- factor(output[["meta_data"]]$sample_id[mm], levels=samplevels)
umapRes$cell_clustering = factor(umapRes$cell_clustering, levels=clusterlevels)

pdf('plot_umaps.pdf',width=10,height=10)
plotUmap(umapRes = umapRes,
         code_clustering=cell_clustering1m,
         color_clusters = clustercolors,
         subtype_markers = output$subtype_markers)
dev.off()

#saveRDS(umapRes, file="backup_umap.rds")


####FUNCTIONAL BOXES####

fmlistplot <- output$functional_markers
exprtbl <- 
  data.frame(fsApply(output$fcs,exprs)[, union(output$ƒsubtype_markers,output$functional_markers)],
             sample_id = output$sample_ids, cluster = output$cell_clustering1m) %>%
  group_by(sample_id, cluster) %>%
  summarize_all(funs(mean))

ggdf2<-melt(exprtbl, id.var=c("cluster","sample_id"))
ggdf2$cluster <- factor(ggdf2$cluster, levels=clusterlevels)
ggdf2$condition<- factor(output$meta_data$condition[match(ggdf2$sample_id,output$meta_data$sample_id)],levels=grouplevels)
ggdf2$sample_id <- factor(ggdf2$sample_id, levels = samplevels)
ggdf2<-ggdf2[ggdf2$cluster %nin% c("UA"),]

#groups for statistics
comps<-list(c("I","IT"),c("LI","LIT"))

library(ggpubr)

pdf("plot_funcmarkers_stats.pdf",width=6,height=7)
for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf2[ggdf2$variable==fmlistplot[i],], aes(x=condition, y=value, fill=condition))+
    geom_boxplot(outlier.size=0, lwd=0.25, color="black")+
    geom_jitter(width=0, size=0.5)+
    facet_wrap(~cluster,ncol=4,scales="free")+
    ylab("MMI")+
    stat_compare_means(method="t.test", label="p.format",method.args = list(var.equal = TRUE),label.x.npc="center",size=1.5,comparisons=comps)+
    ggtitle(fmlistplot[i])+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=9, color="black"),
          axis.title.x = element_blank(),
          axis.line.x = element_line(size=0.25, color="black"),
          axis.line.y = element_line(size=0.25, color="black"),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size=12, color="black"),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=7, color="black"),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key.size = unit(1.5,'lines'),
          legend.text = element_text(size=6),
          legend.key = element_rect(fill="white"))
 
  print(ggp)
}
dev.off()


##end