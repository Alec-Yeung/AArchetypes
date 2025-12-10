rm(list=ls())
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(openxlsx)

###
AA_82_pathways <- read.xlsx('AA_82_pathways.xlsx')
AA_diets <- read_excel("AA_diets.xlsx")

AA_82_pathways_1<- AA_82_pathways %>% 
  mutate(`D+N` = D + N, `E+Q` = E + Q)
AA_82_pathways_1$D <- NULL
AA_82_pathways_1$N <- NULL
AA_82_pathways_1$E <- NULL
AA_82_pathways_1$Q <- NULL

AA_82_pathways_2 <- AA_82_pathways_1 %>% 
  dplyr::select(Cluster, A:`E+Q`)

common_cols <- intersect(colnames(AA_82_pathways_2)[-1], colnames(AA_diets)[-1])

AA_diets <- AA_diets[, c("Diet", common_cols)]
colnames(AA_diets)[1] <- "Cluster"

AA_pathways_diets <- rbind(AA_82_pathways_2, AA_diets)  

### t-SNE
library(Rtsne)

set.seed(11)
tsne_result <- Rtsne(AAs, dims = 2, perplexity = 30)

tsne_plot_data <- as.data.frame(tsne_result$Y)
colnames(tsne_plot_data) <- c("Dim1", "Dim2")
tsne_plot_data$Cluster <- labels

colors <- c("A" = "#f57c6e", "B" = "#fbe79e", "C" = "#f2b56e", "D" = "#b9bddc",
            "Mediterranean" = "#8dd3c7", "Japanese" = "#ff6700", "Vegetarian" = "#bebada",
            "Plant-based" = "#432818", "DASH" = "#80b1d3", "Paleo" = "#fdb462", "Ketogenic" = "#b3de69",
            "Atkins" = "#fccde5", "American" = "#d9d9d9", "USDA" = "#bc80bd")

shapes <- c("A" = 16, "B" = 16, "C" = 16, "D" = 16,
            "Mediterranean" = 2, "Japanese" = 2, "Vegetarian" = 2,
            "Plant-based" = 2, "DASH" = 2, "Paleo" = 2, "Ketogenic" = 2,
            "Atkins" = 2, "American" = 2, "USDA" = 2)

ordered_levels <- c("A", "B", "C", "D", "Mediterranean", "Japanese",
                    "Vegetarian", "Plant-based", "DASH", "Paleo",
                    "Ketogenic", "Atkins", "American", "USDA")

tsne_plot_data$Cluster <- factor(tsne_plot_data$Cluster, levels = ordered_levels)

P <- ggplot(tsne_plot_data, aes(x = Dim1, y = Dim2, color = Cluster, shape = Cluster)) +
  theme_bw()+
  geom_point(size = 3) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  theme_minimal() +
  labs(title = 't-SNE of Pathways+Diets', x = "Dim 1", y = "Dim 2")+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(color = 'black'),
        axis.text.x = element_text(color = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1))
P

ggsave("TSNE of Pathways+Diets2.jpg", width = 6, height = 5, limitsize = FALSE, dpi = 1000)





