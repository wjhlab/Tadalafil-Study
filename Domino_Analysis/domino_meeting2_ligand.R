#combining the ligand from two filtered groups and merging them into full ligand list

tcell1v2 <- read.csv("Tcell_TopDE_1v2unfiltered.csv")
tcell1v2$cell_type <- "Tcell"
tcell4v5 <- read.csv("Tcell_TopDE_4v5unfiltered.csv")
tcell4v5$cell_type <- "Tcell"
gra1v2 <- read.csv("Granulocyte_TopDE_1v2unfiltered.csv")
gra1v2$cell_type <- "Granulocyte"
gra4v5 <- read.csv("Granulocyte_TopDE_4v5unfiltered.csv")
gra4v5$cell_type <- "Granulocyte"
mac1v2 <- read.csv("Macrophage_TopDE_1v2unfiltered.csv")
mac1v2$cell_type <- "Macrophage"
mac4v5 <- read.csv("Macrophage_TopDE_4v5unfiltered.csv")
mac4v5$cell_type <- "Macrophage"

ligand1v2 <- readRDS("signaling1v2.RDS")
list1v2 <- rownames(ligand1v2)
ligand4v5 <- readRDS("signaling4v5.RDS")
list4v5 <- rownames(ligand4v5)

oneVStwo <- rbind(tcell1v2,gra1v2,mac1v2)
fourVSfive <- rbind(tcell4v5,gra4v5,mac4v5)

oneVStwo <- oneVStwo[oneVStwo$X %in% list1v2,]
oneVStwo <- oneVStwo[order(oneVStwo$X),]
oneVStwo$trend <- ifelse(oneVStwo$log2FoldChange>0,'VEH','TAD')
fourVSfive <- fourVSfive[fourVSfive$X %in% list4v5,]
fourVSfive <- fourVSfive[order(fourVSfive$X),]
fourVSfive$trend <- ifelse(fourVSfive$log2FoldChange>0,'MLIST','MLIST.TAD')

write.csv(oneVStwo,"ligandTrend1v2.csv")
write.csv(fourVSfive,"ligandTrend4v5.csv")
