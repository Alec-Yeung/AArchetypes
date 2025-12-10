rm(list=ls())
library(openxlsx)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(ggsignif)

ensembl_entrez_K <- read.xlsx('Evolutionary_rate.xlsx')

proteins <- read.xlsx('Combined_ages_in_pathways.xlsx')

proteins <- proteins %>%
  rename(entrezgene_id = Gene)

ensembl_entrez_K$entrezgene_id <- as.character(ensembl_entrez_K$entrezgene_id)
proteins$entrezgene_id <- as.character(proteins$entrezgene_id)

proteins_K <- proteins %>%
  left_join(ensembl_entrez_K, by = 'entrezgene_id')

proteins_K <- proteins_K %>% 
  mutate(evolutionary_rate = 10 ^ protein_sequence_pecentage_difference)

violin_data <- proteins_K %>%
  dplyr::select(Cluster, evolutionary_rate)

###
P <- ggplot(violin_data, aes(x=as.factor(Cluster), y=evolutionary_rate, 
                             fill=Cluster)) +
  geom_violin(trim=TRUE ,fill="#adb5bd", color="white",
              scale = 'width', width = 0.8, alpha = 1) + 
  geom_boxplot(width=0.1, position=position_dodge(0.9),
               outlier.colour = NA, fill = 'white')+
  theme_bw()+ 
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        axis.title.y=element_text(size = 12, color = 'black'), 
        axis.title.x=element_text(size = 11), 
        axis.text.y=element_text(size=12, color = 'black'),
        axis.text.x=element_text(size=15, color = 'black'),
        legend.title=element_text(size = 15),
        legend.text=element_text(size = 12), 
        panel.border = element_rect(size = 1, colour = 'black'), 
        # axis.line = element_line(colour = "black"), #将x=0轴和y=0轴加粗显示(size=1)
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+
  labs(title = 'Evolutionary rate by AArchetypes', 
       x = '', y = '')+ 
  guides(fill=guide_legend(title = NULL))+
  ylim(1,1.4)+
  stat_compare_means(method = "kruskal.test", label.y = 1.4, size = 5)
P

kruskal.test(evolutionary_rate ~ Cluster, data = violin_data) # Kruskal-Wallis chi-squared = 88.009, df = 3, p-value < 2.2e-16
pchisq(88.009, df = 3, lower.tail = FALSE) # 5.863342e-19

ggsave('Evolutionary rate of clusters of proteins.jpg', P, 
       width = 4, height = 4, limitsize = FALSE, dpi = 1000)

