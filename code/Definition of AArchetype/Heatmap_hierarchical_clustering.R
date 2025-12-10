rm(list=ls())
library(pheatmap)
library(openxlsx)
library(ggplot2)
library(grid)
library(cluster)
library(clValid)
library(dplyr)

data <- read.xlsx('AA_enrichment_score_hsa.xlsx', rowNames = TRUE)

dist_matrix <- dist(data, method = "euclidean")

hc <- hclust(dist_matrix, method = "ward.D2")

k <- 4 ##
clusters <- cutree(hc, k)

sil <- silhouette(clusters, dist_matrix)
plot(sil, border = NA)
mean_silhouette <- mean(sil[, 3]) 

row_annotation <- data.frame(AArchetype = factor(clusters))
rownames(row_annotation) <- rownames(data)

cluster_colors <- list(AArchetype = c("1" = "#f57c6e", "2" = "#f2b56e",
                                      "3" = "#fbe79e", "4" = "#b9bddc")) 
									  
P <- pheatmap(data,
              annotation_row = row_annotation,
              annotation_colors = cluster_colors,
              annotation_legend = FALSE,
              main = "", 
              cluster_rows = hc,  
              color = colorRampPalette(c("#9ecafe", "white", "#fecee6"))(100), 
              legend = TRUE,
              angle_col = 45,
              fontsize = 17,
              fontsize_row = 15,
              fontsize_col = 15,
              treeheight_col = 20,
              cutree_rows = 4, ##
              annotation_names_row = FALSE) 

pdf('AA Enrichment Scores of Human Metabolic Pathways.pdf',
    width = 15, height = 16)
print(P)
grid.text("AA Enrichment Scores of Human Metabolic Pathways", 
          x = 0.3, y = 0.99, gp = gpar(fontface = "plain", fontsize = 20))
grid.text("AArchetype", x = 0.08, y = 0.06, just = "right",
          gp = gpar(fontface = "plain", fontsize = 17), rot = 45)
grid.text("AA enrichment score", x = 0.93, y = 0.89, 
          gp = gpar(fontsize = 22, fontface = "plain"), rot = 90)
grid.text("A", x = 0.06, y = 0.90, 
          gp = gpar(fontsize = 22, fontface = "plain"), rot = 0)
grid.text("B", x = 0.06, y = 0.67, 
          gp = gpar(fontsize = 22, fontface = "plain"), rot = 0)
grid.text("C", x = 0.06, y = 0.44, 
          gp = gpar(fontsize = 22, fontface = "plain"), rot = 0)
grid.text("D", x = 0.06, y = 0.20,
          gp = gpar(fontsize = 22, fontface = "plain"), rot = 0)
dev.off()

###########################################################################
sil_df <- as.data.frame(sil[, 1:3])
rownames(sil_df) <- rownames(data)
# write.xlsx(sil_df, 'Silhouette_coefficients.xlsx', rowNames = TRUE)

sil_df$Pathway <- rownames(sil_df)

sil_df$cluster <- factor(sil_df$cluster, levels = c(1, 2, 3, 4),
                         labels = c("A", "C", "B", "D"))
sil_df$cluster <- factor(sil_df$cluster, levels = c("D", "C", "B", "A"))

sil_df <- sil_df %>%
  arrange(cluster, sil_width)

sil_df$Pathway <- factor(sil_df$Pathway, levels = sil_df$Pathway)

cluster_colors <- c("A" = "#fccde5", "B" = "#bebada",
                    "C" = "#a2d2ff", "D" = "#80a1d4")

P <- ggplot(sil_df, aes(x = reorder(Pathway, cluster), y = sil_width, fill = as.factor(cluster))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = cluster_colors, name = "AArchetype",
                    breaks = c("A", "B", "C", "D"),
                    guide = guide_legend(keywidth = 3, keyheight = 3)) +
  labs(title = "Silhouette coefficients", x = "", y = "") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 17),
    axis.text.x = element_text(size = 17, color = 'black'),
    axis.text.y = element_text(size = 15, color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
P

ggsave('Silhouette coefficients of metabolic pathways.jpg', P, width = 15, height = 16, limitsize = FALSE)




