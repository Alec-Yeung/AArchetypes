rm(list = ls())
library(ggplot2)
library(openxlsx)
library(ggpubr)

data <- read.xlsx('Combined_IDR_in_pathways.xlsx')

violin_data <- data %>%
  dplyr::select(Cluster, IDR_proportion)

P <- ggplot(violin_data, 
            aes(x=Cluster, y=IDR_proportion, fill=Cluster)) + 
  geom_violin(width=0.9, trim=TRUE,fill="#adb5bd",color="white") + 
  geom_boxplot(width=0.08,position=position_dodge(0.9),outlier.colour = NA, fill = 'white')+
  theme_bw()+ 
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.y=element_text(size = 17), 
        axis.text.y=element_text(size=15, color = 'black'),
        axis.text.x=element_text(size=15, color = 'black'),
        legend.title=element_text(size = 15),
        legend.text=element_text(size = 12), 
        panel.border = element_rect(size = 1, color = 'black'), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  labs(title = 'IDR Proportion by AArchetypes', x = '', y = '')+
  guides(fill=guide_legend(title = NULL))+
  coord_cartesian(ylim = c(0, 0.4))+
  stat_compare_means(comparisons = list(c("A", "B"),
                                        c("C", "D")),
                     method = "wilcox.test", label = "p.signif", size = 5, bracket.size = 0.7,
                     tip.length = 0.01,
                     label.y = c(.3, .3)
  )
P     

ggsave('IDR Proportion by AArchetypes.jpg', P, width = 4, height = 4, limitsize = FALSE, dpi = 1000)


