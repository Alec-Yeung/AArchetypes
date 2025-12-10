rm(list=ls())
library(readxl)
library(dplyr)
library(tidyr)

file_paths <- list.files(path = "Protein_intensity/", pattern = "\\.xlsx$", full.names = TRUE)

protein_data <- list()

for (file in file_paths) {
  protein_id <- tools::file_path_sans_ext(basename(file))
  
  df <- read_excel(file)
  
  # 只保留需要的3列 (ID、组织名称、表达量)
  df_clean <- df %>%
    select(UNIQUE_IDENTIFIER, TISSUE_NAME, NORMALIZED_INTENSITY)

  # 将数据存储到列表中
  protein_data[[protein_id]] <- df_clean
}

# 合并所有数据到一个数据框
combined_data <- bind_rows(protein_data)

# 转换为宽格式数据，构建矩阵
expression_matrix <- combined_data %>%
  pivot_wider(
    names_from = TISSUE_NAME, 
    values_from = NORMALIZED_INTENSITY,
    values_fill = NA # 如果某蛋白在某组织没有数据，填充为 NA
  )

# save(expression_matrix, file = 'Matrix_Protein-Tissue.Rdata')
# write.xlsx(as.data.frame(expression_matrix), file = "Matrix_Protein-Tissue.xlsx")

##################### Protein Expression -> Gene Expression ####################
library(biomaRt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

uniprot_id <- expression_matrix$UNIQUE_IDENTIFIER

uniprot_2_entrez <- getBM(attributes = c('uniprotsptrembl', 'entrezgene_id'),
                   filters = 'uniprotsptrembl',
                   values = uniprot_id,
                   mart = ensembl)
colnames(uniprot_2_entrez)[1] <- 'UNIQUE_IDENTIFIER'

expression_matrix1 <- expression_matrix %>% 
  inner_join(uniprot_2_entrez, by = "UNIQUE_IDENTIFIER")

expression_matrix1 <- expression_matrix1 %>%
  mutate(across(`adrenal gland`:`aqueous humour`, as.numeric))  

expression_matrix2 <- expression_matrix1 %>%
  group_by(entrezgene_id) %>%                    
  summarise(across(`adrenal gland`:`aqueous humour`, ~ sum(.x, na.rm = TRUE)))        

# save(expression_matrix2, file = 'Matrix_Gene-Tissue.Rdata')
# write.xlsx(as.data.frame(expression_matrix2), file = "Matrix_Gene-Tissue.xlsx")

##################### Tissue-Specific Correlation ##############################
Combined_costs_of_proteins <- read.xlsx('Combined_costs_of_proteins.xlsx')

cost_data <- Combined_costs_of_proteins %>% dplyr::select(Gene, Protein_cost, Energy_cost)
colnames(cost_data)[1] <- 'entrezgene_id'
cost_data <- cost_data %>%
  distinct(entrezgene_id, .keep_all = TRUE)

expression_matrix3 <- expression_matrix2 %>% 
  inner_join(cost_data, by = 'entrezgene_id')

expression_matrix3 <- expression_matrix3 %>%
  mutate(across(`adrenal gland`:Energy_cost, as.numeric))  




corr_df_P <- data.frame(
  Tissue = character(),  
  Correlation = numeric(),  
  P_Value = numeric()  
)

for (col_name in names(expression_matrix3)[2:79]) {  
  result <- cor.test(expression_matrix3[[col_name]], expression_matrix3$Protein_cost, method = "spearman")  
  
  corr_df_P <- rbind(
    corr_df_P,
    data.frame(
      Tissue = col_name,
      Correlation = result$estimate,
      P_Value = result$p.value
    )
  )
}
rownames(corr_df_P) <- NULL

# write.xlsx(corr_df_P, file = "Corr_TissueSpecific_ProteinCost.xlsx")


corr_df_E <- data.frame(
  Tissue = character(),  
  Correlation = numeric(),  
  P_Value = numeric()  
)

for (col_name in names(expression_matrix3)[2:79]) {  
  result <- cor.test(expression_matrix3[[col_name]], expression_matrix3$Energy_cost, method = "spearman")  
  
  corr_df_E <- rbind(
    corr_df_E,
    data.frame(
      Tissue = col_name,
      Correlation = result$estimate,
      P_Value = result$p.value
    )
  )
}
rownames(corr_df_E) <- NULL

# write.xlsx(corr_df_E, file = "Corr_TissueSpecific_EnergyCost.xlsx")

############################## Tissue * Cluster ################################
cluster_gene <- Combined_costs_of_proteins %>%
  dplyr::select(Cluster, Gene)
colnames(cluster_gene)[2] <- 'entrezgene_id'

cluster_gene$entrezgene_id <- as.character(cluster_gene$entrezgene_id)
expression_matrix3$entrezgene_id <- as.character(expression_matrix3$entrezgene_id)

expression_matrix4 <- cluster_gene %>% 
  left_join(expression_matrix3, by = 'entrezgene_id')
expression_matrix4 <- expression_matrix4[!apply(expression_matrix4[, 3:80], 1, function(row) all(is.na(row))), ]
expression_matrix4$Protein_cost <- NULL
expression_matrix4$Energy_cost <- NULL
expression_matrix4$entrezgene_id <- NULL

expression_matrix5 <- expression_matrix4 %>%
  group_by(Cluster) %>%                          
  summarise(across(`adrenal gland`:`aqueous humour`, sum)) 
expression_matrix5 <- t(expression_matrix5)
colnames(expression_matrix5) <- expression_matrix5[1,]
expression_matrix5 <- expression_matrix5[-1,]
expression_matrix5 <- as.data.frame(expression_matrix5)
# write.xlsx(expression_matrix5, 'Tissue_Cluster.xlsx', rowNames=TRUE)
expression_matrix5 <- type.convert(expression_matrix5, as.is = TRUE)


library(reldist)
gini_scores <- apply(expression_matrix5, 2, gini)
print(gini_scores)


cv_scores <- apply(expression_matrix5, 2, function(x) sd(x) / mean(x))
print(cv_scores)


result1 <- cor.test(expression_matrix5$A, expression_matrix5$D, method = "spearman") 
correlation <- result1$estimate; correlation
p_value <- result1$p.value; p_value







