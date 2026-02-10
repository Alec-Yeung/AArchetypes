rm(list=ls())
library(randomForest)
library(dplyr)
library(caret)
library(pROC)
library(openxlsx)

data <- read.xlsx("Combined_genes_in_pathways.xlsx")  

data_RF <- data %>%
  dplyr::select(-Pathway_ID, -Pathway, -Gene)

data_RF <- data_RF %>%
  dplyr::mutate(Cluster = ifelse(Cluster == "D", 1, 0)) ##

set.seed(42)
train_index <- createDataPartition(data_RF$Cluster, p = 0.7, list = FALSE)
train_data <- data_RF[train_index,]
test_data <- data_RF[-train_index,]

rf_model <- randomForest(as.factor(Cluster) ~ ., data = train_data,
                         ntree = 500, mtry = 4, importance = TRUE)
print(rf_model)

rf_predictions <- predict(rf_model, newdata = test_data)

conf_matrix <- confusionMatrix(rf_predictions, as.factor(test_data$Cluster))
print(conf_matrix)

rf_prob_predictions <- predict(rf_model, newdata = test_data,
                               type = "prob")[,2]

pdf('RF_D-NOT D.pdf', width = 10, height = 10, pointsize = 20) ##
roc_curve <- roc(test_data$Cluster, rf_prob_predictions, levels = c("0", "1"),
                 plot=TRUE, auc.polygon=FALSE,
                 print.auc=TRUE, show.thres=TRUE,
                 legacy.axes=TRUE, print.thres=FALSE,
                 xlab='FPR', ylab='TPR', main='') 
title(main = "D and NOT D", font.main = 1, line = 2.2) ##
dev.off()



