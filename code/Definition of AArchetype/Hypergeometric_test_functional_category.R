rm(list=ls())
library(openxlsx)
library(ggplot2)

df <- read.xlsx('Pathways_functional_mapping.xlsx', rowNames = TRUE) 

clusters <- unique(df$cluster) 
categories <- unique(df$Functional_category) 

p_matrix <- matrix(NA, nrow = length(clusters), ncol = length(categories))
rownames(p_matrix) <- clusters
colnames(p_matrix) <- categories

total_paths <- nrow(df) 

for (cluster in clusters) {
  n_A <- sum(df$cluster == cluster)  
  
  for (category in categories) {
    n_category_total <- sum(df$Functional_category == category)  
    n_category_A <- sum(df$cluster == cluster & df$Functional_category == category)  
    
    p_value <- phyper(n_category_A - 1, n_category_total, total_paths - n_category_total, n_A, lower.tail = FALSE)
    
    p_matrix[cluster, category] <- p_value
  }
}

print(p_matrix)

df_p_matrix <- as.data.frame(p_matrix)
write.xlsx(df_p_matrix, 'p_hyper.xlsx', rowNames = TRUE) 


