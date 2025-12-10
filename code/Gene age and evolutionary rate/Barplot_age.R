rm(list=ls())
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)

data <- read.xlsx('Combined_ages_in_pathways.xlsx')
data <- data[!is.na(data$Age), ]

data1 <- subset(data, Age %in% c('454.6','361.2','324.5','220.2','176.1','104.7'))

data2 <- data1 %>% group_by(Cluster, Age) %>% tally(name = 'Count')

data3 <- data2 %>% group_by(Cluster) %>% 
  mutate(Proportion = Count/sum(Count))
data3$Age <- factor(data3$Age)

# Chi-square test of independence / Chi-square goodness-of-fit test
count_cb <- xtabs(Count ~ Cluster + Age,
                  data = data3[data3$Cluster %in% c("C", "B"), ])
chisq_test_cb <- chisq.test(count_cb)
p_cb <- formatC(chisq_test_cb$p.value, format = "e", digits = 1)

count_cd <- xtabs(Count ~ Cluster + Age,
                  data = data3[data3$Cluster %in% c("C", "D"), ])
chisq_test_cd <- chisq.test(count_cd)
p_cd <- formatC(chisq_test_cd$p.value, format = "e", digits = 1)

count_ca <- xtabs(Count ~ Cluster + Age,
                  data = data3[data3$Cluster %in% c("C", "A"), ])
chisq_test_ca <- chisq.test(count_ca)
p_ca <- formatC(chisq_test_ca$p.value, format = "e", digits = 1)



P <- ggplot(data3, aes(x = Cluster, y = Proportion, fill=Age))+ 
  geom_bar(stat = 'identity', width = 0.7)+
  ggtitle('Age distribution of AArchetypes')+
  xlab('')+
  ylab('Proportion')+
  labs(fill = 'Age (mya)')+
  scale_y_continuous(limits = c(0,1.15), breaks=seq(0, 1, by=0.2))+
  scale_fill_manual(values = c("#e76253","#fbb75f",
                               "#fde7b9","#85c6d8",
                               "#4783a5","#1e476d"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5),
        legend.title = element_text(size = 12, color = 'black'),
        panel.border = element_rect(colour = 'black', size = 2), 
        axis.text.x = element_text(size = 15, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 11, color = 'black'),
        axis.title.x = element_text(size = 12, color = 'black'),
        axis.title.y = element_text(size = 12, color = 'black'))+
  geom_signif(comparisons = list(c("A", "B"), c("C", "D")),
              tip_length = 0.02,
              size = 1,
              textsize = 4.5,
              annotations = c('****', '****'),
              y_position = c(1.00, 1.00)
              )
P

ggsave('Age distribution of AArchetypes.jpg',P,
       width = 5, height = 6, limitsize = FALSE, dpi = 1000)

