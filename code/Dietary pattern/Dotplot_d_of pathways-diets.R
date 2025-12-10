rm(list=ls())
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(tidyr)
library(dplyr)

###
data <- read.xlsx('Distance_pathways_diets.xlsx')
data$Gene <- NULL

data1 <- data %>%
  pivot_longer(
    cols = -Cluster,        
    names_to = "Diet",      
    values_to = "d"         
  )

medians <- data1 %>%
  group_by(Diet, Cluster) %>%
  summarise(median = median(d, na.rm = TRUE)) %>%
  ungroup()

medians$Cluster <- factor(medians$Cluster, levels = c("A", "B", "C", "D"))
medians$Diet <- factor(medians$Diet)

###
P <- ggplot(medians, aes(x = Cluster, y = Diet,
                        color = Cluster, size = median)) +
  geom_point(alpha = 0.8, shape = 16) +
  scale_size_continuous(range = c(14, 18), breaks = c(14, 16, 18)) +
  scale_color_manual(values = c("A" = "#f57c6e", 
                                "B" = "#fbe79e", 
                                "C" = "#f2b56e", 
                                "D" = "#b9bddc")) +
  theme_minimal() +
  labs(title = "Distance between Diets and Pathways in AArchetypes",
       x = "",
       y = "",
       color = "AArchetype",
       size = "Median") +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.text.x = element_text(size = 14, color = 'black', hjust = 1),
    axis.text.y = element_text(size = 14, color = 'black'),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
P




P <- ggplot(medians, aes(x = Cluster, y = Diet,
                         color = Cluster, size = median)) +
  geom_point(alpha = 0.8, shape = 16) +
  scale_size_continuous(range = c(5, 15),  
                        breaks = c(14, 15, 16, 17, 18)) +
  scale_color_manual(values = c("A" = "#f57c6e", 
                                "B" = "#fbe79e", 
                                "C" = "#f2b56e", 
                                "D" = "#b9bddc")) +
  theme_minimal() +
  labs(title = "Distance between Diets and Pathways in AArchetypes",
       x = "",
       y = "",
       color = "AArchetype",
       size = "Median") +  
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.text.x = element_text(size = 15, color = 'black'),
    axis.text.y = element_text(size = 15, color = 'black'),
    legend.title = element_text(size = 13, colour = 'black'),
    legend.text = element_text(size = 12, colour = 'black'),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  geom_text(aes(label = round(median, 1)), size = 2, color = "black")+
  guides(color = guide_legend(override.aes = list(size = 5)))

P




