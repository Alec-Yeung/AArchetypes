rm(list=ls())
library(dplyr)
library(reldist)
library(openxlsx)

load('Matrix_Gene-Tissue.Rdata')
sum(is.na(expression_matrix2$entrezgene_id))
expression_matrix2 <- expression_matrix2 %>% filter(!is.na(entrezgene_id))
expression_matrix2 <- as.data.frame(expression_matrix2)
rownames(expression_matrix2) <- expression_matrix2$entrezgene_id
expression_matrix2$entrezgene_id <- NULL

cluster_genes <- read.xlsx('Combined_genes_in_pathways.xlsx')
cluster_genes <- cluster_genes %>%
  dplyr::select(Cluster, Gene)
cluster_genes$Gene <- as.character(cluster_genes$Gene)

################################# Gini #########################################
gini_scores <- apply(expression_matrix2, 1, gini)
names(gini_scores) <- rownames(expression_matrix2)
head(gini_scores)

df_gini <- as.data.frame(gini_scores)
df_gini$Gene <- rownames(df_gini)
rownames(df_gini) <- NULL

merge_df_gini <- cluster_genes %>%
  inner_join(df_gini, by = 'Gene')

############################ Violinplot: Gini ##################################
library(ggplot2)
library(ggpubr)

violin_data_Gini <- merge_df_gini %>%
  dplyr::select(Cluster, gini_scores)

P <- ggplot(violin_data_Gini, 
            aes(x=Cluster, y=gini_scores, fill=Cluster)) + 
  geom_violin(width=0.8, trim=TRUE,fill="#adb5bd",color="white") + 
  geom_boxplot(width=0.1,position=position_dodge(0.9),outlier.colour = NA, fill = 'white')+
  theme_bw()+ 
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.y=element_text(size = 17), 
        axis.text.y=element_text(size=15, color = 'black'),
        axis.text.x=element_text(size=15, color = 'black'),
        legend.title=element_text(size = 15),
        legend.text=element_text(size = 12), 
        panel.border = element_rect(size = 1, colour = 'black'), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  labs(title = 'Gini scores by AArchetypes', x = '', y = '')+
  guides(fill=guide_legend(title = NULL))+
  stat_compare_means(comparisons = list(c("A", "B"),
                                        c("C", "D")
                                        ),
                     method = "wilcox.test", label = "p.signif", size = 4, bracket.size = 0.7, vjust = 0.4,
                     tip.length = 0.02,
                     label.y = c(1, 1)
  )
P     

ggsave('Gini scores by AArchetypes.jpg', P, width = 4, height = 4, limitsize = FALSE, dpi = 1000)

