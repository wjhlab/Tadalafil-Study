#Nicole Gross
#R version 4.2.3 (2023-03-15)
#Platform: aarch64-apple-darwin20 (64-bit)
#Running under: macOS 14.0

rm(list = ls())
library(reshape2)
library(randomcoloR)
library(pals)
library(ggplot2)
library(Hmisc)
library(edgeR)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(corrplot)
library(readxl)
library(ggridges)
library(ComplexHeatmap)
library(limma)
library(FlowSOM)
library(ggpubr)
library(matrixStats)
library(dplyr)

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
  md$mouse <- factor(md$mouse)
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
  sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))
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
    exprs(x) <- expr
    x
  })
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
  som <- ReadInput(fcs, transform = FALSE, scale = FALSE) %>% BuildSOM(colsToUse = subtype_markers)
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
                                           cluster_merging = NULL, 
                                            subtype_markers,
                                            clusterMergeFile=NULL,
                                            fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  color_clusters <-c(stepped2(20),stepped3(20),rev(cubicl(20)))
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
  color_clusters=clustercolors
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
                               samplevels = samplevels,
                               subtype_markers = output$subtype_markers,
                               color_conditions = clustercolors,
                               shape_conditions = c(1:13),
                               fileName = 'diagnostics.pdf', 
                               tit = '', 
                               fun = mean)
{
  pdf(file = fileName)
  
  ## Spot check - number of cells per sample
  cell_table <- table(sample_ids)
  ggdf <- data.frame(sample_id = names(cell_table), 
                     cell_counts = as.numeric(cell_table))
  ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
  mm <- match(ggdf$sample_id, md$sample_id)
  ggdf$condition <- md$condition[mm]
  print(ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = condition)) + 
          geom_bar(stat = 'identity') + 
          geom_text(aes(label = cell_counts), hjust = 0.5, vjust = -0.5, size = 2.5) + 
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  
          scale_fill_manual(values = color_conditions, drop = FALSE) + 
          scale_x_discrete(drop = FALSE))
  
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
                    clusterMergeFile='~/Desktop/ViralHCC/ViralHCC_merging.xlsx',
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
  ####umapindex generation####
  #umap ncells = table of sample ids with how many to downsample to by default col = id, row = ncells
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

plotUmap <- function(umapRes,
                     seed=1234,
                     neighbors=10,
                     midpoint,
                     color_clusters,
                     code_clustering,subtype_markers=NULL)
{require(umap);require(ggplot2);require(viridis);require(ggrepel)
  custom.settings = umap.defaults
  custom.settings$seed = seed
  custom.settings$n.neighbors = neighbors
  ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = cell_clustering)) +
    geom_point(size = 1) +
    #geom_text(aes(label=cell_clustering), size=1) +
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    
    scale_color_manual(values = color_clusters, name="CLUSTERS") +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  
  print(ggp)
  
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

###Skip lines 447-457 if loading output rds
#read and cluster
output <- returnfcs(metaDataFile = paste0(workd,"/Config","/metadata.xlsx"),
                    panelDataFile = paste0(workd,"/Config","/panel.xlsx"),
                    dataDirectory = paste0(workd,'/Data'))

output[8:10] <- clusterfcs(numclusters=40)

names(output)[8:10] <- c('code_clustering','cell_clustering','metaclusters')

#assign clusters
clusterMergeFile = paste0(workd,"/Config",'/merge.xlsx')
cluster_merging <- read_xlsx(clusterMergeFile)

#clusterlevels=c(1:40) #assigning cluster levels before annotation
clusterlevels = c("B",
                  "cDC",
                  "cDC2",
                  "Endoth",
                  "Gran_I",
                  "Gran_II",
                  "Mac_I",
                  "Mac_II",
                  "Mac_III",
                  "Mac_IV",
                  "Monocyte_I",
                  "Monocyte_II",
                  "Monocyte_III",
                  "NI",
                  "NK",
                  "Stroma",
                  "Tc_I",
                  "Tc_II",
                  "Th_I",
                  "Th_II",
                  "Th_III",
                  "Tumor_I",
                  "Tumor_II",
                  "UA")

samplevels <- c("mus1", "mus2", "mus3", "mus4", "mus5",
                "mus6", "mus7", "mus8", "mus9", "mus10",
                "mus11", "mus12", "mus13", "mus14", "mus15",
                "mus16", "mus17", "mus18", "mus19", "mus20",
                "mus21", "mus22", "mus23", "mus24", "mus25",
                "mus26", "mus27", "mus28", "mus29", "mus30",
                "mus31", "mus32", "mus33", "mus34", "mus35",
                "mus36", "mus37", "mus38", "mus39", "mus40",
                "mus41", "mus42")

condlevels=c("VEH",
             "TAD")

batchlevels=c("1", "2", "3", "4", "5")

grouplevels = c("grp1", "grp2", "grp3", "grp4", "grp5", "grp6")

combotxlevels=c("ISO", "PD1", "CTLA4")

txlevels=c("VEH_ISO","TAD_ISO","VEH_PD1","TAD_PD1","VEH_CTLA4","TAD_CTLA4")

clustercolors <- as.character(c(cols25(n=25),alphabet(n=19)))

###Skip lines 509-511 if loading output rds
mm1 <- match(output$cell_clustering, cluster_merging$original_cluster)
cell_clustering1m <- cluster_merging$new_cluster[mm1]
output$cell_clustering1m <- cell_clustering1m

#metacluster heatmap
plot_clustering_heatmap_wrapper2(fcs=output$fcs,
                                 color_clusters = clustercolors,
                                 cell_clustering = factor(output$cell_clustering1m, levels=clusterlevels), 
                                 subtype_markers=output$subtype_markers,
                                 fileName = 'clusteringheatmap_merged.pdf');dev.off()

#save output
#saveRDS(output, file = "backup_output.rds")

expr <- fsApply(output$fcs, exprs) #create expression matrix
rng <- colQuantiles(expr, probs = c(0.01, 0.99))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1])) #scaling 0-1
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1 # need this expr object later 


####DIFFERENTIAL PLOTS####

#set up count and prop matrices for immune cells only
counts_table_imm <- table(output$cell_clustering1m[output$cell_clustering1m %nin% c("Endoth","NI","Stroma","Tumor_I","Tumor_II","UA")], 
                          output$sample_ids[output$cell_clustering1m %nin% c("Endoth","NI","Stroma","Tumor_I","Tumor_II","UA")])
props_table_imm <- t(t(counts_table_imm) / colSums(counts_table_imm)) * 100
counts_imm <- as.data.frame.matrix(counts_table_imm)
props_imm <- as.data.frame.matrix(props_table_imm)

write.csv(counts_imm, file='countsimm.csv')
write.csv(props_imm, file='propsimm.csv')

#set up the data frame for proportional plotting of immune cells
ggdf <- melt(data.frame(cluster = rownames(props_imm),props_imm, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdf$cluster <- factor(ggdf$cluster, levels=clusterlevels)
ggdf$condition <- factor(output$meta_data$condition[match(ggdf$sample_id,output$meta_data$sample_id)], levels=condlevels)
ggdf$batch <- output$meta_data$batch[match(ggdf$sample_id,output$meta_data$sample_id)]
ggdf$combotx <- factor(output$meta_data$combotx[match(ggdf$sample_id,output$meta_data$sample_id)], levels=combotxlevels)
ggdf$group <- factor(output$meta_data$group[match(ggdf$sample_id,output$meta_data$sample_id)], levels=grouplevels)
ggdf$tx <- factor(output$meta_data$tx[match(ggdf$sample_id,output$meta_data$sample_id)], levels=txlevels)

grp1to4 <- ggdf[ggdf$group %in% c("grp1", "grp2", "grp3", "grp4"),]
grp1to4$group <- factor(grp1to4$group, levels=c("grp1", "grp2", "grp3", "grp4"))

#heatmap of select markers for TAM clusters
expr01 <-expr01[,output$subtype_markers[output$subtype_markers!="PDPN"]]
expr01_mean <- data.frame(expr01, cell_clustering = output$cell_clustering1m, check.names = FALSE) %>%
  group_by(cell_clustering) %>% summarize_all(funs(mean))
expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
rownames(expr_heat) <- expr01_mean$cell_clustering
expr_heat_TAM <- expr_heat[c("Mac_I","Mac_II","Mac_III","Mac_IV"),
                           c("F480","CD19","CD62L","IAIE","CD11C","PDL1", "Ly6C","CD206")]
pdf("plot_heatmap_TAM.pdf",width=4,height=4)
Heatmap(expr_heat_TAM, name="scaled expr",
        width = ncol(expr_heat_TAM)*unit(5, "mm"), 
        height = nrow(expr_heat_TAM)*unit(5, "mm"))
dev.off()

#groups for statistics 
comps<-list(c("grp1","grp2"),c("grp3","grp4"))

#plot box plots: immune cells only for group 1-4 
boxplot <- ggplot(grp1to4, aes(x=group, y=proportion, fill=group))+
  geom_boxplot(outlier.size=0, lwd=0.25, color="black")+
  geom_jitter(width=0, size=0.5)+
  facet_wrap(~cluster,ncol=8,scales="free")+
  ylab("% of Immune Cells")+
  stat_compare_means(method="t.test", label="p.format",method.args = list(var.equal = TRUE),label.x.npc="center", size=1.5,aes(group = group),comparisons= comps)+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=7, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white")
  )
pdf("plot_abundanceimm_box_grp1to4_stats.pdf", width=9, height=5);boxplot;dev.off()


## Umap

#to load umap rds
umapRes<-readRDS('backup_umap.rds')

###skip lines 606-617 if loading UMAP rds
umapRes <- do_umap(fcs=output$fcs,subtype_markers = output$subtype_markers,
                   sample_ids = output$sample_ids,cell_clustering = output$cell_clustering, metadata=output$metadata,
                   clusterMergeFile=clusterMergeFile,
                   seed = 1234, ncells=500,sample_subset=NULL)

mm <- match(as.character(umapRes$sample_id), as.character(output[["meta_data"]]$sample_id))
umapRes$condition <- factor(output[["meta_data"]]$condition[mm], levels=condlevels)
umapRes$group <- factor(output[["meta_data"]]$group[mm], levels=grouplevels)
umapRes$batch <- output[["meta_data"]]$batch[mm]
umapRes$sample_id <- factor(output[["meta_data"]]$sample_id[mm], levels=samplevels)
umapRes$cell_clustering = factor(umapRes$cell_clustering, levels=clusterlevels)
umapRes<-umapRes[umapRes$cell_clustering!="UA",]

dev.off()
pdf('plot_umaps.pdf',width=10,height=10)
plotUmap(umapRes = umapRes,
         code_clustering=cell_clustering1m,
         color_clusters = clustercolors,
         subtype_markers = output$subtype_markers)
dev.off()

#saveRDS(umapRes, file="backup_umap.rds")


## Functional markers
markerlist = output$functional_markers

exprtbl <-
  data.frame(fsApply(output$fcs,exprs)[, union(output$subtype_markers,output$functional_markers)],
             sample_id = output$sample_ids, cluster = output$cell_clustering1m) %>%
  group_by(sample_id, cluster) %>%
  summarize_all(funs(mean))
ggdf2<-melt(exprtbl, id.var=c("cluster","sample_id"))
ggdf2$cluster <- factor(ggdf2$cluster, levels=clusterlevels)
ggdf2$condition <- factor(output$meta_data$condition[match(ggdf2$sample_id,output$meta_data$sample_id)], levels=condlevels)
ggdf2$combotx <- factor(output$meta_data$combotx[match(ggdf2$sample_id,output$meta_data$sample_id)], levels=combotxlevels)
ggdf2$sample_id <- factor(ggdf2$sample_id, levels = samplevels)
ggdf2$group <- factor(output$meta_data$group[match(ggdf2$sample_id,output$meta_data$sample_id)], levels=grouplevels)

grp1to4<-ggdf2[ggdf2$group %in% c("grp1","grp2","grp3","grp4"),]
grp1to4$group<-factor(grp1to4$group, levels=grouplevels)

#Box plots
pdf("plot_box_funcmarkers_grp1to4_stats.pdf",width=9,height=5)
for(i in 1:length(markerlist)){
  ggp <- ggplot(grp1to4[grp1to4$variable==markerlist[i],], aes(x=group, y=value, fill=group))+
    geom_boxplot(outlier.size=0, lwd=0.25, color="black")+
    geom_jitter(width=0, size=0.5)+
    facet_wrap(~cluster,ncol=8,scales="free")+
    ylab("MMI")+
    ggtitle(markerlist[i])+
    stat_compare_means(method="t.test", label="p.format",method.args = list(var.equal = TRUE),label.x.npc="center", size=1.5,aes(group = group),comparisons= comps)+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=10, color="black"),
          axis.title.x = element_blank(),
          axis.line.x = element_line(size=0.25, color="black"),
          axis.line.y = element_line(size=0.25, color="black"),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size=12, color="black"),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=7, color="black"),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key.size = unit(2,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white")
    )
  print(ggp)
}
dev.off()

####end