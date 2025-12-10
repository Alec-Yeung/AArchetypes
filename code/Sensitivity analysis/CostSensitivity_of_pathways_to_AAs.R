rm(list=ls())
library(dplyr)
library(openxlsx)

### AA_pathways
AA_pathways <- read.xlsx('AA_82_pathways.xlsx')

AA_pathways[, 4:23] <- AA_pathways[, 4:23]/100

AA_protein_cost <- c(0.300433333,
                     8.661673333,
                     1.13525,
                     24.44967667,
                     0,
                     22.49612667,
                     0,
                     0,
                     0,
                     0,
                     0,
                     28.79824,
                     18.45821667,
                     11.55789333,
                     25.98524333,
                     110.8193,
                     0,
                     0,
                     0,
                     3.9915
)

AA_energy_cost <- c(0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    5,
                    0,
                    5,
                    15,
                    0,
                    0,
                    0,
                    0,
                    0
)

AA_pathways$Protein_cost <- as.matrix(AA_pathways[, 4:23]) %*% AA_protein_cost
AA_pathways$Energy_cost <- as.matrix(AA_pathways[, 4:23]) %*% AA_energy_cost

AA_data <- as.matrix(AA_pathways[4:23])

dx <- 0.001

### Cost sensitivity to AAs - Protein cost
s_matrix1 <- matrix(NA, nrow = 82, ncol = 20)
colnames(s_matrix1) <- colnames(AA_data)
rownames(s_matrix1) <- AA_pathways$Pathway

for (i in 1:nrow(s_matrix1)) {
  AA_data_i <- AA_data[i,] 
  
  for (j in 1:ncol(s_matrix1)){
    AA_data_i[j] <- AA_data_i[j] + dx 
    
    c_after <- (AA_data_i %*% AA_protein_cost) / (1+dx)
    c_after <- as.numeric(c_after)
    c_before <- AA_pathways$Protein_cost[i] 
    
    s <- (c_after - c_before) / dx
    
    s_matrix1[i, j] <- s
    
    AA_data_i[j] <- AA_data_i[j] - dx 
  }
}

# write.xlsx(s_matrix1, file = 'CostSensitivity_of_pathways_to_AAs_p.xlsx', rowNames = TRUE)

### Cost sensitivity to AAs - Energy cost
s_matrix2 <- matrix(NA, nrow = 82, ncol = 20)
colnames(s_matrix2) <- colnames(AA_data)
rownames(s_matrix2) <- AA_pathways$Pathway

for (i in 1:nrow(s_matrix2)) {
  AA_data_i <- AA_data[i,] 
  
  for (j in 1:ncol(s_matrix2)){
    AA_data_i[j] <- AA_data_i[j] + dx 
    
    c_after <- (AA_data_i %*% AA_energy_cost) / (1+dx)
    c_after <- as.numeric(c_after)
    c_before <- AA_pathways$Energy_cost[i] 
    
    s <- (c_after - c_before) / dx
    
    s_matrix2[i, j] <- s
    
    AA_data_i[j] <- AA_data_i[j] - dx 
  }
}

# write.xlsx(s_matrix2, file = 'CostSensitivity_of_pathways_to_AAs_e.xlsx', rowNames = TRUE)

### Convert s_matrix to -1/1 matrix
# 使用 apply 逐行处理
s_matrix_transform <- t(apply(s_matrix, 1, function(row) {
  # 初始化为 0
  res <- rep(0, length(row))
  
  # 标记最小值位置为 -1
  res[row == min(row, na.rm = TRUE)] <- -1
  
  # 标记最大值位置为 1
  res[row == max(row, na.rm = TRUE)] <- 1
  
  return(res)
}))

# 设置行名和列名
rownames(s_matrix_transform) <- rownames(s_matrix)
colnames(s_matrix_transform) <- colnames(s_matrix)

# write.xlsx(s_matrix_transform, file = 'CostSensitivity_of_pathways_to_AAs_transform.xlsx', rowNames = TRUE)





