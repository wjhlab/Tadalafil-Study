library(corrplot)
library(ggcorrplot)
library(GGally)
library(ggplot2)
library(tidyverse)
library(dplyr)
setwd("E:/WJHLab/revised_domino/tirangle_heatmap")
vehtad_hm <- read.csv("vehtad_hm.csv")
vehtad_hm <- vehtad_hm[,!colnames(vehtad_hm)=='X']
vehtad_hm <- vehtad_hm[vehtad_hm$receptor %in% c("Csf1r",
                                             "Adora2a",
                                             "Ccr5",
                                             "Ccr2",
                                             "Il6ra",
                                             "Il4ra",
                                             "Mrc1",
                                             "Notch1",
                                             "Notch2",
                                             "Il10rb",
                                             "Ltbr"),]
colnames(vehtad_hm)[1] <- 'X'
data_matrix <- as.matrix(vehtad_hm)
data_matrix1 <- data_matrix[,c(2:9)]

# Get dimensions of the data matrix
n_rows <- nrow(data_matrix1)
n_cols <- ncol(data_matrix1)

# Generate shifted coordinates for polygons
polygonsblue <- do.call(rbind,
                     apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                       row_idx = xy[1]
                       col_idx = xy[2]
                       cell_value = data_matrix1[row_idx, col_idx]
                       a <- row_idx - 1
                       b <- col_idx - 1
                       if (cell_value %in% c(2, 3)) {
                         data.frame(
                           x = c(0, 0, 1) + b,#upper triangle coordinates
                           y = c(0, 1, 1) + a,#upper triangle coordinates
                           group = paste(a, b, sep = "-"),
                           value = cell_value
                         )
                       } else {
                         NULL
                       }
                     }))

polygonsyellow <- do.call(rbind,
                     apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                       row_idx = xy[1]
                       col_idx = xy[2]
                       cell_value = data_matrix1[row_idx, col_idx]
                       a <- row_idx - 1
                       b <- col_idx - 1
                       if (cell_value %in% c(1, 3)) {
                         data.frame(
                           x = c(0, 1, 1) + b,#lower triangle coordinates
                           y = c(0, 0, 1) + a,#lower triangle coordinates
                           group = paste(a, b, sep = "-"),
                           value = cell_value
                         )
                       } else {
                         NULL
                       }
                     }))

polygonswhiteup <- do.call(rbind,
                     apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                       row_idx = xy[1]
                       col_idx = xy[2]
                       cell_value = data_matrix1[row_idx, col_idx]
                       a <- row_idx - 1
                       b <- col_idx - 1
                       if (cell_value %in% c(0,1)) {
                         data.frame(
                           x = c(0, 0, 1) + b,#upper triangle coordinates
                           y = c(0, 1, 1) + a,#upper triangle coordinates
                           group = paste(a, b, sep = "-"),
                           value = cell_value
                         )
                       } else {
                         NULL
                       }
                     }))

polygonswhitelower <- do.call(rbind,
                           apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                             row_idx = xy[1]
                             col_idx = xy[2]
                             cell_value = data_matrix1[row_idx, col_idx]
                             a <- row_idx - 1
                             b <- col_idx - 1
                             if (cell_value %in% c(0,2)) {
                               data.frame(
                                 x = c(0, 1, 1) + b,#lower triangle coordinates
                                 y = c(0, 0, 1) + a,#lower triangle coordinates
                                 group = paste(a, b, sep = "-"),
                                 value = cell_value
                               )
                             } else {
                               NULL
                             }
                           }))


pdf("test.pdf",height = 35, width = 10)
plot1 = ggplot() +
  geom_polygon(data = polygonsblue, aes(x, y, group = group), color = "black",fill = "blue") +
  geom_polygon(data = polygonsyellow, aes(x, y, group = group), color = "black",fill = "orange") +
  geom_polygon(data = polygonswhiteup, aes(x, y, group = group), color = "black",fill = "white") +
  geom_polygon(data = polygonswhitelower, aes(x, y, group = group), color = "black",fill = "white") +
  theme_void() +
  coord_fixed(ratio = 1) +
  scale_x_continuous(expand = c(0, 0), breaks = 1:n_cols, labels = colnames(data_matrix1)) +
  scale_y_continuous(expand = c(0, 0), breaks = 1:n_rows, labels = vehtad_hm$X) +
  labs(x = "Column Title", y = "Row Title",fill = "Custom Legend Title") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = -1, hjust = 0),
        axis.text.y = element_text(angle = 0, vjust = 2, hjust = 0),
        plot.margin = margin(1, 1, 1, 1, "cm"))
print(plot1)
dev.off()


mlistmlisttad_hm=read.csv("mlistmlisttad_hm.csv")
mlistmlisttad_hm <- mlistmlisttad_hm[,!colnames(mlistmlisttad_hm)=='X']
mlistmlisttad_hm <- mlistmlisttad_hm[mlistmlisttad_hm$receptor %in% c("Notch1","Notch2","Cxcr3","Il12rb2"),]
colnames(mlistmlisttad_hm)[1] <- 'X'
data_matrix2 <- as.matrix(mlistmlisttad_hm)
data_matrix3 <- data_matrix2[,c(2:9)]

# Get dimensions of the data matrix
n_rows <- nrow(data_matrix3)
n_cols <- ncol(data_matrix3)

# Generate shifted coordinates for polygons
polygonsblue1 <- do.call(rbind,
                        apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                          row_idx = xy[1]
                          col_idx = xy[2]
                          cell_value = data_matrix3[row_idx, col_idx]
                          a <- row_idx - 1
                          b <- col_idx - 1
                          if (cell_value %in% c(1, 3)) {
                            data.frame(
                              x = c(0, 0, 1) + b,#upper triangle coordinates
                              y = c(0, 1, 1) + a,#upper triangle coordinates
                              group = paste(a, b, sep = "-"),
                              value = cell_value
                            )
                          } else {
                            NULL
                          }
                        }))

polygonsyellow1 <- do.call(rbind,
                          apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                            row_idx = xy[1]
                            col_idx = xy[2]
                            cell_value = data_matrix3[row_idx, col_idx]
                            a <- row_idx - 1
                            b <- col_idx - 1
                            if (cell_value %in% c(2, 3)) {
                              data.frame(
                                x = c(0, 1, 1) + b,#lower triangle coordinates
                                y = c(0, 0, 1) + a,#lower triangle coordinates
                                group = paste(a, b, sep = "-"),
                                value = cell_value
                              )
                            } else {
                              NULL
                            }
                          }))

polygonswhiteup1 <- do.call(rbind,
                           apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                             row_idx = xy[1]
                             col_idx = xy[2]
                             cell_value = data_matrix3[row_idx, col_idx]
                             a <- row_idx - 1
                             b <- col_idx - 1
                             if (cell_value %in% c(0,2)) {
                               data.frame(
                                 x = c(0, 0, 1) + b,#upper triangle coordinates
                                 y = c(0, 1, 1) + a,#upper triangle coordinates
                                 group = paste(a, b, sep = "-"),
                                 value = cell_value
                               )
                             } else {
                               NULL
                             }
                           }))

polygonswhitelower1 <- do.call(rbind,
                              apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                                row_idx = xy[1]
                                col_idx = xy[2]
                                cell_value = data_matrix3[row_idx, col_idx]
                                a <- row_idx - 1
                                b <- col_idx - 1
                                if (cell_value %in% c(0,1)) {
                                  data.frame(
                                    x = c(0, 1, 1) + b,#lower triangle coordinates
                                    y = c(0, 0, 1) + a,#lower triangle coordinates
                                    group = paste(a, b, sep = "-"),
                                    value = cell_value
                                  )
                                } else {
                                  NULL
                                }
                              }))

pdf("test1.pdf",height = 35, width = 10)
plot2 = ggplot() +
  geom_polygon(data = polygonsblue1, aes(x, y, group = group), color = "black",fill = "blue") +
  geom_polygon(data = polygonsyellow1, aes(x, y, group = group), color = "black",fill = "orange") +
  geom_polygon(data = polygonswhiteup1, aes(x, y, group = group), color = "black",fill = "white") +
  geom_polygon(data = polygonswhitelower1, aes(x, y, group = group), color = "black",fill = "white") +
  theme_void() +
  coord_fixed(ratio = 1) +
  scale_x_continuous(expand = c(0, 0), breaks = 1:n_cols, labels = colnames(data_matrix3)) +
  scale_y_continuous(expand = c(0, 0), breaks = 1:n_rows, labels = mlistmlisttad_hm$X) +
  labs(x = "Column Title", y = "Row Title") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = -1, hjust = 0),
        axis.text.y = element_text(angle = 0, vjust = 2, hjust = 0),
        plot.margin = margin(1, 1, 1, 1, "cm"))+
  guides(fill = guide_legend(title = "Legend Title"))
print(plot2)
dev.off()

library(gridExtra)
library(grid)

pdf("full_plotss.pdf",height = 30, width = 10)
combined_plot <- grid.arrange(
  plot1 + labs(subtitle = "VEH vs TAD") +
    theme(plot.margin = margin(10, 10, 10, 10, "pt")),
  plot2 + labs(subtitle = "MLIST vs MLIST.TAD") +
    theme(plot.margin = margin(10, 10, 10, 10, "pt")),
  ncol = 2,  top = textGrob("Heatmap for receptors by different comparisons", gp = gpar(fontsize = 20)
))
dev.off()




###########################################keep rows only have 1
vehtad_hm=read.csv("vehtad_hm.csv")
vehtad_hm <- vehtad_hm[,!colnames(vehtad_hm)=='X']
colnames(vehtad_hm)[1] <- 'X'
filtered_data <- vehtad_hm %>%
  filter(rowSums(. == 1) > 0)
filtered_data1 <- filtered_data[,c(2:9)]
data_matrix_filtered <- as.matrix(filtered_data1)
# Get dimensions of the data matrix
n_rows <- nrow(data_matrix_filtered)
n_cols <- ncol(data_matrix_filtered)

# Generate shifted coordinates for polygons
polygonsblue <- do.call(rbind,
                        apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                          row_idx = xy[1]
                          col_idx = xy[2]
                          cell_value = data_matrix_filtered[row_idx, col_idx]
                          a <- row_idx - 1
                          b <- col_idx - 1
                          if (cell_value %in% c(2, 3)) {
                            data.frame(
                              x = c(0, 0, 1) + b,#upper triangle coordinates
                              y = c(0, 1, 1) + a,#upper triangle coordinates
                              group = paste(a, b, sep = "-"),
                              value = cell_value
                            )
                          } else {
                            NULL
                          }
                        }))

polygonsyellow <- do.call(rbind,
                          apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                            row_idx = xy[1]
                            col_idx = xy[2]
                            cell_value = data_matrix_filtered[row_idx, col_idx]
                            a <- row_idx - 1
                            b <- col_idx - 1
                            if (cell_value %in% c(1, 3)) {
                              data.frame(
                                x = c(0, 1, 1) + b,#lower triangle coordinates
                                y = c(0, 0, 1) + a,#lower triangle coordinates
                                group = paste(a, b, sep = "-"),
                                value = cell_value
                              )
                            } else {
                              NULL
                            }
                          }))

polygonswhiteup <- do.call(rbind,
                           apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                             row_idx = xy[1]
                             col_idx = xy[2]
                             cell_value = data_matrix_filtered[row_idx, col_idx]
                             a <- row_idx - 1
                             b <- col_idx - 1
                             if (cell_value %in% c(0,1)) {
                               data.frame(
                                 x = c(0, 0, 1) + b,#upper triangle coordinates
                                 y = c(0, 1, 1) + a,#upper triangle coordinates
                                 group = paste(a, b, sep = "-"),
                                 value = cell_value
                               )
                             } else {
                               NULL
                             }
                           }))

polygonswhitelower <- do.call(rbind,
                              apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                                row_idx = xy[1]
                                col_idx = xy[2]
                                cell_value = data_matrix_filtered[row_idx, col_idx]
                                a <- row_idx - 1
                                b <- col_idx - 1
                                if (cell_value %in% c(0,2)) {
                                  data.frame(
                                    x = c(0, 1, 1) + b,#lower triangle coordinates
                                    y = c(0, 0, 1) + a,#lower triangle coordinates
                                    group = paste(a, b, sep = "-"),
                                    value = cell_value
                                  )
                                } else {
                                  NULL
                                }
                              }))


pdf("test.pdf",height = 35, width = 10)
plot1 = ggplot() +
  geom_polygon(data = polygonsblue, aes(x, y, group = group), color = "black",fill = "blue") +
  geom_polygon(data = polygonsyellow, aes(x, y, group = group), color = "black",fill = "orange") +
  geom_polygon(data = polygonswhiteup, aes(x, y, group = group), color = "black",fill = "white") +
  geom_polygon(data = polygonswhitelower, aes(x, y, group = group), color = "black",fill = "white") +
  theme_void() +
  coord_fixed(ratio = 1) +
  scale_x_continuous(expand = c(0, 0), breaks = 1:n_cols, labels = colnames(data_matrix_filtered)) +
  scale_y_continuous(expand = c(0, 0), breaks = 1:n_rows, labels = filtered_data$X) +
  labs(x = "Column Title", y = "Row Title",fill = "Custom Legend Title") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = -1, hjust = 0),
        axis.text.y = element_text(angle = 0, vjust = 2, hjust = 0),
        plot.margin = margin(1, 1, 1, 1, "cm"))
print(plot1)
dev.off()


mlistmlisttad_hm=read.csv("mlistmlisttad_hm.csv")
mlistmlisttad_hm <- mlistmlisttad_hm[,!colnames(mlistmlisttad_hm)=='X']
colnames(mlistmlisttad_hm)[1] <- 'X'
filtered_data <- mlistmlisttad_hm %>%
  filter(rowSums(. == 2) > 0)
filtered_data1 <- filtered_data[,c(2:9)]
data_matrix_filtered <- as.matrix(filtered_data1)

# Get dimensions of the data matrix
n_rows <- nrow(data_matrix_filtered)
n_cols <- ncol(data_matrix_filtered)

# Generate shifted coordinates for polygons
polygonsblue1 <- do.call(rbind,
                         apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                           row_idx = xy[1]
                           col_idx = xy[2]
                           cell_value = data_matrix_filtered[row_idx, col_idx]
                           a <- row_idx - 1
                           b <- col_idx - 1
                           if (cell_value %in% c(1, 3)) {
                             data.frame(
                               x = c(0, 0, 1) + b,#upper triangle coordinates
                               y = c(0, 1, 1) + a,#upper triangle coordinates
                               group = paste(a, b, sep = "-"),
                               value = cell_value
                             )
                           } else {
                             NULL
                           }
                         }))

polygonsyellow1 <- do.call(rbind,
                           apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                             row_idx = xy[1]
                             col_idx = xy[2]
                             cell_value = data_matrix_filtered[row_idx, col_idx]
                             a <- row_idx - 1
                             b <- col_idx - 1
                             if (cell_value %in% c(2, 3)) {
                               data.frame(
                                 x = c(0, 1, 1) + b,#lower triangle coordinates
                                 y = c(0, 0, 1) + a,#lower triangle coordinates
                                 group = paste(a, b, sep = "-"),
                                 value = cell_value
                               )
                             } else {
                               NULL
                             }
                           }))

polygonswhiteup1 <- do.call(rbind,
                            apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                              row_idx = xy[1]
                              col_idx = xy[2]
                              cell_value = data_matrix_filtered[row_idx, col_idx]
                              a <- row_idx - 1
                              b <- col_idx - 1
                              if (cell_value %in% c(0,2)) {
                                data.frame(
                                  x = c(0, 0, 1) + b,#upper triangle coordinates
                                  y = c(0, 1, 1) + a,#upper triangle coordinates
                                  group = paste(a, b, sep = "-"),
                                  value = cell_value
                                )
                              } else {
                                NULL
                              }
                            }))

polygonswhitelower1 <- do.call(rbind,
                               apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                                 row_idx = xy[1]
                                 col_idx = xy[2]
                                 cell_value = data_matrix_filtered[row_idx, col_idx]
                                 a <- row_idx - 1
                                 b <- col_idx - 1
                                 if (cell_value %in% c(0,1)) {
                                   data.frame(
                                     x = c(0, 1, 1) + b,#lower triangle coordinates
                                     y = c(0, 0, 1) + a,#lower triangle coordinates
                                     group = paste(a, b, sep = "-"),
                                     value = cell_value
                                   )
                                 } else {
                                   NULL
                                 }
                               }))

pdf("test1.pdf",height = 35, width = 10)
plot2 = ggplot() +
  geom_polygon(data = polygonsblue1, aes(x, y, group = group), color = "black",fill = "blue") +
  geom_polygon(data = polygonsyellow1, aes(x, y, group = group), color = "black",fill = "orange") +
  geom_polygon(data = polygonswhiteup1, aes(x, y, group = group), color = "black",fill = "white") +
  geom_polygon(data = polygonswhitelower1, aes(x, y, group = group), color = "black",fill = "white") +
  theme_void() +
  coord_fixed(ratio = 1) +
  scale_x_continuous(expand = c(0, 0), breaks = 1:n_cols, labels = colnames(data_matrix_filtered)) +
  scale_y_continuous(expand = c(0, 0), breaks = 1:n_rows, labels = filtered_data$X) +
  labs(x = "Column Title", y = "Row Title") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = -1, hjust = 0),
        axis.text.y = element_text(angle = 0, vjust = 2, hjust = 0),
        plot.margin = margin(1, 1, 1, 1, "cm"))+
  guides(fill = guide_legend(title = "Legend Title"))
print(plot2)
dev.off()

library(gridExtra)
library(grid)

pdf("full_plot_filtered.pdf",height = 30, width = 10)
combined_plot <- grid.arrange(
  plot1 + labs(subtitle = "VEH vs TAD") +
    theme(plot.margin = margin(10, 10, 10, 10, "pt")),
  plot2 + labs(subtitle = "MLIST vs MLIST.TAD") +
    theme(plot.margin = margin(10, 10, 10, 10, "pt")),
  ncol = 2,  top = textGrob("Heatmap for receptors by different comparisons", gp = gpar(fontsize = 20)
  ))
dev.off()






###########################################keep single blue and single orange
vehtad_hm=read.csv("vehtad_hm.csv")
vehtad_hm <- vehtad_hm[,!colnames(vehtad_hm)=='X']
colnames(vehtad_hm)[1] <- 'X'
filtered_data <- vehtad_hm[rowSums(vehtad_hm == 3) == 0, ]
filtered_data <- filtered_data[rowSums(filtered_data == 0) != 8, ]
filtered_data1 <- filtered_data[,c(2:9)]
data_matrix_filtered <- as.matrix(filtered_data1)
# Get dimensions of the data matrix
n_rows <- nrow(data_matrix_filtered)
n_cols <- ncol(data_matrix_filtered)

# Generate shifted coordinates for polygons
polygonsblue <- do.call(rbind,
                        apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                          row_idx = xy[1]
                          col_idx = xy[2]
                          cell_value = data_matrix_filtered[row_idx, col_idx]
                          a <- row_idx - 1
                          b <- col_idx - 1
                          if (cell_value %in% c(2, 3)) {
                            data.frame(
                              x = c(0, 0, 1) + b,#upper triangle coordinates
                              y = c(0, 1, 1) + a,#upper triangle coordinates
                              group = paste(a, b, sep = "-"),
                              value = cell_value
                            )
                          } else {
                            NULL
                          }
                        }))

polygonsyellow <- do.call(rbind,
                          apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                            row_idx = xy[1]
                            col_idx = xy[2]
                            cell_value = data_matrix_filtered[row_idx, col_idx]
                            a <- row_idx - 1
                            b <- col_idx - 1
                            if (cell_value %in% c(1, 3)) {
                              data.frame(
                                x = c(0, 1, 1) + b,#lower triangle coordinates
                                y = c(0, 0, 1) + a,#lower triangle coordinates
                                group = paste(a, b, sep = "-"),
                                value = cell_value
                              )
                            } else {
                              NULL
                            }
                          }))

polygonswhiteup <- do.call(rbind,
                           apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                             row_idx = xy[1]
                             col_idx = xy[2]
                             cell_value = data_matrix_filtered[row_idx, col_idx]
                             a <- row_idx - 1
                             b <- col_idx - 1
                             if (cell_value %in% c(0,1)) {
                               data.frame(
                                 x = c(0, 0, 1) + b,#upper triangle coordinates
                                 y = c(0, 1, 1) + a,#upper triangle coordinates
                                 group = paste(a, b, sep = "-"),
                                 value = cell_value
                               )
                             } else {
                               NULL
                             }
                           }))

polygonswhitelower <- do.call(rbind,
                              apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                                row_idx = xy[1]
                                col_idx = xy[2]
                                cell_value = data_matrix_filtered[row_idx, col_idx]
                                a <- row_idx - 1
                                b <- col_idx - 1
                                if (cell_value %in% c(0,2)) {
                                  data.frame(
                                    x = c(0, 1, 1) + b,#lower triangle coordinates
                                    y = c(0, 0, 1) + a,#lower triangle coordinates
                                    group = paste(a, b, sep = "-"),
                                    value = cell_value
                                  )
                                } else {
                                  NULL
                                }
                              }))


pdf("test.pdf",height = 35, width = 10)
plot1 = ggplot() +
  geom_polygon(data = polygonsblue, aes(x, y, group = group), color = "black",fill = "blue") +
  geom_polygon(data = polygonsyellow, aes(x, y, group = group), color = "black",fill = "orange") +
  geom_polygon(data = polygonswhiteup, aes(x, y, group = group), color = "black",fill = "white") +
  geom_polygon(data = polygonswhitelower, aes(x, y, group = group), color = "black",fill = "white") +
  theme_void() +
  coord_fixed(ratio = 1) +
  scale_x_continuous(expand = c(0, 0), breaks = 1:n_cols, labels = colnames(data_matrix_filtered)) +
  scale_y_continuous(expand = c(0, 0), breaks = 1:n_rows, labels = filtered_data$X) +
  labs(x = "Column Title", y = "Row Title",fill = "Custom Legend Title") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = -1, hjust = 0),
        axis.text.y = element_text(angle = 0, vjust = 2, hjust = 0),
        plot.margin = margin(1, 1, 1, 1, "cm"))
print(plot1)
dev.off()


mlistmlisttad_hm=read.csv("mlistmlisttad_hm.csv")
mlistmlisttad_hm <- mlistmlisttad_hm[,!colnames(mlistmlisttad_hm)=='X']
colnames(mlistmlisttad_hm)[1] <- 'X'
filtered_data <- mlistmlisttad_hm[rowSums(mlistmlisttad_hm == 3) == 0, ]
filtered_data <- filtered_data[rowSums(filtered_data == 0) != 8, ]
filtered_data1 <- filtered_data[,c(2:9)]
data_matrix_filtered <- as.matrix(filtered_data1)

# Get dimensions of the data matrix
n_rows <- nrow(data_matrix_filtered)
n_cols <- ncol(data_matrix_filtered)

# Generate shifted coordinates for polygons
polygonsblue1 <- do.call(rbind,
                         apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                           row_idx = xy[1]
                           col_idx = xy[2]
                           cell_value = data_matrix_filtered[row_idx, col_idx]
                           a <- row_idx - 1
                           b <- col_idx - 1
                           if (cell_value %in% c(1, 3)) {
                             data.frame(
                               x = c(0, 0, 1) + b,#upper triangle coordinates
                               y = c(0, 1, 1) + a,#upper triangle coordinates
                               group = paste(a, b, sep = "-"),
                               value = cell_value
                             )
                           } else {
                             NULL
                           }
                         }))

polygonsyellow1 <- do.call(rbind,
                           apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                             row_idx = xy[1]
                             col_idx = xy[2]
                             cell_value = data_matrix_filtered[row_idx, col_idx]
                             a <- row_idx - 1
                             b <- col_idx - 1
                             if (cell_value %in% c(2, 3)) {
                               data.frame(
                                 x = c(0, 1, 1) + b,#lower triangle coordinates
                                 y = c(0, 0, 1) + a,#lower triangle coordinates
                                 group = paste(a, b, sep = "-"),
                                 value = cell_value
                               )
                             } else {
                               NULL
                             }
                           }))

polygonswhiteup1 <- do.call(rbind,
                            apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                              row_idx = xy[1]
                              col_idx = xy[2]
                              cell_value = data_matrix_filtered[row_idx, col_idx]
                              a <- row_idx - 1
                              b <- col_idx - 1
                              if (cell_value %in% c(0,2)) {
                                data.frame(
                                  x = c(0, 0, 1) + b,#upper triangle coordinates
                                  y = c(0, 1, 1) + a,#upper triangle coordinates
                                  group = paste(a, b, sep = "-"),
                                  value = cell_value
                                )
                              } else {
                                NULL
                              }
                            }))

polygonswhitelower1 <- do.call(rbind,
                               apply(expand.grid(1:n_rows, 1:n_cols), 1, function(xy) {
                                 row_idx = xy[1]
                                 col_idx = xy[2]
                                 cell_value = data_matrix_filtered[row_idx, col_idx]
                                 a <- row_idx - 1
                                 b <- col_idx - 1
                                 if (cell_value %in% c(0,1)) {
                                   data.frame(
                                     x = c(0, 1, 1) + b,#lower triangle coordinates
                                     y = c(0, 0, 1) + a,#lower triangle coordinates
                                     group = paste(a, b, sep = "-"),
                                     value = cell_value
                                   )
                                 } else {
                                   NULL
                                 }
                               }))

pdf("test1.pdf",height = 15, width = 10)
plot2 = ggplot() +
  geom_polygon(data = polygonsblue1, aes(x, y, group = group), color = "black",fill = "blue") +
  geom_polygon(data = polygonsyellow1, aes(x, y, group = group), color = "black",fill = "orange") +
  geom_polygon(data = polygonswhiteup1, aes(x, y, group = group), color = "black",fill = "white") +
  geom_polygon(data = polygonswhitelower1, aes(x, y, group = group), color = "black",fill = "white") +
  theme_void() +
  coord_fixed(ratio = 1) +
  scale_x_continuous(expand = c(0, 0), breaks = 1:n_cols, labels = colnames(data_matrix_filtered)) +
  scale_y_continuous(expand = c(0, 0), breaks = 1:n_rows, labels = filtered_data$X) +
  labs(x = "Column Title", y = "Row Title") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = -1, hjust = 0),
        axis.text.y = element_text(angle = 0, vjust = 2, hjust = 0),
        plot.margin = margin(1, 1, 1, 1, "cm"))+
  guides(fill = guide_legend(title = "Legend Title"))
print(plot2)
dev.off()

library(gridExtra)
library(grid)

pdf("full_plot_1blue1orange.pdf",height = 30, width = 10)
combined_plot <- grid.arrange(
  plot1 + labs(subtitle = "VEH vs TAD") +
    theme(plot.margin = margin(10, 10, 10, 10, "pt")),
  plot2 + labs(subtitle = "MLIST vs MLIST.TAD") +
    theme(plot.margin = margin(10, 10, 10, 10, "pt")),
  ncol = 2,  top = textGrob("Heatmap for receptors by different comparisons", gp = gpar(fontsize = 20)
  ))
dev.off()