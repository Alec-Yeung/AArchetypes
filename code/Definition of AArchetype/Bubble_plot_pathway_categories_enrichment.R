rm(list=ls())
library(ggplot2)
library(tidyr)
library(dplyr)
library(openxlsx)

p_value_data <- data.frame(
  Carbohydrate = c(0.002446249, 0.811868243, 0.653043306, 0.99053093),
  AA = c(0.642566574, 0.000131773, 0.927615309, 1),
  Global = c(0.291052718, 0.09886084, 0.735130148, 1),
  Lipid = c(0.890515888, 0.998906622, 0.333683627, 0.021215527),
  Cofactors_vitamins = c(1, 0.293995512, 0.735130148, 0.617803893),
  Nucleotide = c(1, 0.122252334, 1, 1),
  Energy = c(1, 0.735478771, 0.381887986, 0.750451671),
  Xenobiotics = c(1, 1, 0.00248419, 1),
  Glycan = c(1, 1, 1, 6.6141e-05),
  Terpenoids_polyketides = c(1, 1, 1, 0.365853659)
)

p_value_data_long <- p_value_data %>%
  mutate(AArchetypes = c('A', 'B', 'C', 'D')) %>%
  pivot_longer(cols = -AArchetypes, names_to = 'Functional_category', values_to = 'p_value')

p_value_data_long1 <- subset(p_value_data_long, p_value < 0.05)

functional_category_data <- read.xlsx('Pathways_functional_mapping.xlsx')
functional_category_data <- functional_category_data %>%
  dplyr::select(Cluster, Functional_category, category_count)
functional_category_data <- unique(functional_category_data)

colnames(functional_category_data)[1] <- 'AArchetypes'

merge_df1 <- p_value_data_long1 %>%
  left_join(functional_category_data, by = c('Functional_category', 'AArchetypes'))

merge_df1 <- merge_df1 %>%
  mutate(`-log10 p-value` = -log10(p_value))

P <- ggplot(merge_df1, aes(x = AArchetypes, y = Functional_category, size = category_count, color = `-log10 p-value`)) +
  geom_point(alpha = ifelse(merge_df1$category_count == 0, 0, 1),  
             shape = ifelse(merge_df1$category_count == 0, NA, 16)) + 
  scale_size_continuous(range = c(3, 15), breaks = c(4,8,12)) +  
  scale_color_gradientn(colors = c("#ebce81", "#e78451", "#aa4744"), 
                        values = c(0,0.25,0.5,0.75,1),  
                        breaks = c(1,2,3,4) 
                        )+ 
  theme_minimal() +  
  labs(title = "Major Categories of Metabolic Pathways \n Enriched in AArchetypes",
       x = "",
       y = "Functional Category",
       color = "-log10 (p-value)",
       size = "Count") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title.y=element_text(size = 13),
        axis.text.x = element_text(size = 12, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'),
        panel.grid.major = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1))
P

ggsave('Major Categories of Metabolic Pathways Enriched in AArchetypes.jpg', P,
       width = 5, height = 5, limitsize = FALSE, dpi = 1000)












