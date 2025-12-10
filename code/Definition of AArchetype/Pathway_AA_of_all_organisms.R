rm(list = ls())
library(openxlsx)
library(ggplot2)
library(dplyr)
library(cowplot)

hsa <- read.xlsx("D:\\Desktop\\AA_composition\\AA_pathways_filterd.xlsx") 
sce <- read.xlsx('sce/AA_pathways_filterd.xlsx') 
ptr <- read.xlsx('ptr/AA_pathways_filterd.xlsx') 
mtu <- read.xlsx('mtu/AA_pathways_filterd.xlsx') 
mmu <- read.xlsx('mmu/AA_pathways_filterd.xlsx') 
eco <- read.xlsx('eco/AA_pathways_filterd.xlsx') 
dme <- read.xlsx('dme/AA_pathways_filterd.xlsx') 
cel <- read.xlsx('cel/AA_pathways_filterd.xlsx') 
bsu <- read.xlsx('bsu/AA_pathways_filterd.xlsx') 
ath <- read.xlsx('ath/AA_pathways_filterd.xlsx') 

df_combined <- bind_rows(hsa, sce, ptr, mtu, mmu, eco, dme, cel, bsu, ath)

df_combined$Organism <- substr(df_combined[[1]], 1, 3)

# write.xlsx(df_combined, 'Pathway_AA_of_organisms.xlsx')

############################### MDS ###################################
### Data preparation
category <- df_combined$Organism
mds_data <- df_combined %>%
  dplyr::select(A:Y)

distance_matrix <- dist(mds_data, method = "euclidean")

mds_result <- cmdscale(distance_matrix, k = 2, eig = TRUE)

mds_coords <- as.data.frame(mds_result$points)
colnames(mds_coords) <- c("Dim1", "Dim2")
head(mds_coords)

mds_plot_data <- cbind(mds_coords, Organism = category)

mds_plot_data <- mds_plot_data %>%
  mutate(Group = case_when(
    Organism %in% c("ath", "cel", "dme", "hsa", "mmu", "ptr") ~ "Multicellular organisms",
    Organism %in% c("bsu", "eco", "mtu", "sce") ~ "Unicellular organisms"
  ))

colors <- c("ath" = "#8dd3c7", "hsa" = "#ff6700", "cel" = "#432818",
            "dme" = "#bebada", "eco" = "#80b1d3", "bsu" = "#fdb462", "mmu" = "#b3de69",
            "mtu" = "#fccde5", "ptr" = "#d9d9d9", "sce" = "#bc80bd")
shapes <- c("mtu" = 8, "eco" = 8, "bsu" = 8, "sce" = 8,
            "ath" = 16, "cel" = 16, "dme" = 16,
            "hsa" = 16, "mmu" = 16, "ptr" = 16)

ordered_levels <- c("ath", "cel", "dme", "hsa", "mmu", "ptr",
                    "bsu", "eco", "mtu", "sce")

mds_plot_data$Organism <- factor(mds_plot_data$Organism, levels = ordered_levels)

P <- ggplot(mds_plot_data, aes(x = Dim1, 
                               y = Dim2,
                               color = Organism,
                               shape = Organism)) +
  theme_bw()+
  geom_point(size = 3) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  theme_minimal() +
  labs(title = "Classical MDS of AA Composition in Organisms' Pathways",
       x = paste('Dim 1'),
       y = paste('Dim 2'),
       color = 'Organism')+
  theme(plot.title = element_text(hjust = 0.3, size = 12),
        legend.title = element_blank(),
        # legend.position = "none",
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 12, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 2))
P

ggsave("Classical MDS of AA Composition in Organisms' Pathways1.jpg", P, width = 5, height = 5, limitsize = FALSE, dpi = 1000)



