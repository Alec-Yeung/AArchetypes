rm(list=ls())
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(tidyr)

ProIntensity_data <- read.xlsx('Combined_ProIntensity_in_pathways.xlsx')

violinplot_data_pi <- ProIntensity_data %>% dplyr::select(Cluster, Protein_intensity)
violinplot_data_pi <- violinplot_data_pi[!is.na(violinplot_data_pi$Protein_intensity), ]

P <- ggplot(violinplot_data_pi,
            aes(x=Cluster, y=Protein_intensity, fill=Cluster)) + 
  geom_violin(width=1.0, trim=TRUE,fill="#adb5bd",color="white") + 
  geom_boxplot(width=0.07,position=position_dodge(0.9),outlier.colour = NA, fill = 'white')+
  theme_bw()+ 
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        axis.title.y=element_text(size = 15), 
        axis.text.y=element_text(size=15, color = 'black'),
        axis.text.x=element_text(size=15, color = 'black'),
        legend.title=element_text(size = 15),
        legend.text=element_text(size = 12), 
        panel.border = element_rect(size = 1, color = 'black'), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+     
  labs(title = 'Protein Intensity by AArchetypes', x = '', y = '')+
  guides(fill=guide_legend(title = NULL))+
  coord_cartesian(ylim = c(0, 3000))+
  stat_compare_means(comparisons = list(c("A", "B"),
                                        c("C", "D")
                                        ),
                     method = "wilcox.test", label = "p.signif", size = 5, bracket.size = 0.7,
                     tip.length = 0.01,
                     label.y = c(2000, 2000)
  )
P      

ggsave('Protein Intensity by AArchetypes.jpg', P, width = 4, height = 4, limitsize = FALSE, dpi = 1000)






