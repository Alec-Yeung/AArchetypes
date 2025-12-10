rm(list=ls())
library(openxlsx)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(pROC)

proteome_LLPS_2_AA <- read.xlsx('proteome_LLPS_2_AA, Human.xlsx')

###
set.seed(111)  
n <- nrow(proteome_LLPS_2_AA)
train_idx <- sample(1:n, size = 0.7*n)

train <- proteome_LLPS_2_AA[train_idx, ]
test  <- proteome_LLPS_2_AA[-train_idx, ]

### 拟合逻辑回归
model <- glm(Is_LLPS ~ LLPS_AA, 
             data = train, 
             family = binomial)

summary(model)

### 画拟合曲线
P <- ggplot(train, aes(x = LLPS_AA, y = Is_LLPS)) +
  geom_point(alpha = 0.3, color = "black") +   # 原始散点
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"),
              se = TRUE, color = "#cf5c4f", fill = "#e29d95") +
  labs(x = "LLPS_AA",
       y = "Probability of LLPS",
       title = "") +
  theme_bw()+
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.y=element_text(size = 15), 
        axis.title.x=element_text(size = 15), 
        axis.text.y=element_text(size=14, color = 'black'),
        axis.text.x=element_text(size=14, color = 'black'),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())

P

ggsave('Logistic curve of LLPS.pdf',
       P, device = cairo_pdf,, width = 7, height = 7,
       limitsize = FALSE, dpi = 1000)


###
test$pred_prob <- predict(model, newdata=test, type="response")

roc_curve <- roc(test$Is_LLPS, test$pred_prob, plot=TRUE,
                 xlab='FPR', ylab='TPR', main='ROC (Test set)')
auc(roc_curve)

###
pdf('ROC curve for LLPS prediction1.pdf', width = 10, height = 10, pointsize = 20) 
roc_curve <- roc(response = test$Is_LLPS, 
                 predictor = test$pred_prob, 
                 levels = c("0", "1"),
                 plot=TRUE, auc.polygon=FALSE,
                 print.auc=TRUE, show.thres=TRUE,
                 legacy.axes=TRUE, print.thres=FALSE,
                 xlab='FPR', ylab='TPR', main='ROC curve for LLPS prediction') 
dev.off()










