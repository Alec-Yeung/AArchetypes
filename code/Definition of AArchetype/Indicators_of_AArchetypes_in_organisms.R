rm(list = ls())
library(openxlsx)
library(ggplot2)
library(dplyr)

data <- read.xlsx('Organism_features.xlsx')

data1 <- data %>%
  mutate(`-log10 (Max(min enrichment p-value))` = -log10(`Max(min.enrichment.p-value)`),
         `-log10 (Average silhouette coefficient)` = -log10(Average_silhouette_coefficient))

P <- ggplot(data1, aes(x = Average_silhouette_coefficient, y = `-log10 (Max(min enrichment p-value))`, 
                      size = Ratio_of_significant_part, color = Is_multicellular_organism)) +
  geom_point(alpha = 1, shape = 16) +
  geom_text(aes(label = Organism,
                hjust = case_when(
                  Organism == "eco" ~ -0.4,
                  Organism == "bsu" ~ 1.3,
                  Organism == "hsa" ~ 1.6,
                  Organism == "dme" ~ -0.2,
                  Organism == "sce" ~ 0.8,
                  Organism == "mmu" ~ -0.2,
                  Organism == "ath" ~ 0.3,
                  TRUE ~ 0.5), 
                vjust = case_when(
                  Organism == "eco" ~ 0.2,
                  Organism == "bsu" ~ 0.2,
                  Organism == "hsa" ~ 1.8,
                  Organism == "dme" ~ 0.5,
                  Organism == "mmu" ~ -0,
                  Organism == "ath" ~ -2.2,
                  Organism == "sce" ~ 1.5,
                  TRUE ~ -1)),
  size = 7, color = 'black', fontface = 'italic') +  
  scale_size_continuous(name = "Ratio of significant part",
                        range = c(2,10), breaks = c(0.25,0.35,0.45)) +  
  scale_color_manual(values = c("Unicellular organism" = "#ffc8c8", "Multicellular organism" = "#cf5149")) +  
  labs(title = "Indicators of AArchetypes in organisms",
       x = "Average silhouette coefficient",
       y = "-log10 (Max(min enrichment p-value))",
       color = "Organism type") +
  theme_minimal()+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x=element_text(size = 14),
        axis.title.y=element_text(size = 14),
        axis.text.x = element_text(size = 13, color = 'black'),
        axis.text.y = element_text(size = 13, color = 'black'),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black", size = 1)
P

ggsave('Indicators of AArchetypes in organisms.pdf',P, width = 8, height = 7,
       bg = "transparent", limitsize = FALSE)










