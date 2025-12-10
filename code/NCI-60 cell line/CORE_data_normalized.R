rm(list = ls())
library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(gg.gap)

CORE_data <- read.xlsx('CORE_data_AAs.xlsx')

CORE_data_long <- CORE_data %>%
  pivot_longer(cols = - Metabolite, 
               names_to = "CellLine", 
               values_to = "Flux")

### except for Gln
CORE_except_gln <- CORE_data_long %>%
  filter(Metabolite != 'Gln')

CORE_except_gln <- CORE_except_gln %>%
  mutate(Metabolite = fct_reorder(Metabolite, Flux, .fun = median))

###
CORE_except_gln1 <- CORE_except_gln %>%
  filter(!Metabolite %in% c('Leu','Lys','Ile','Val','Thr','Phe','Met','Trp'))

P <- ggplot(CORE_except_gln1, 
            aes(x=Metabolite, y=Flux, fill=Metabolite)) + 
  # geom_violin(width=1, trim=TRUE,fill="#adb5bd",color="white") + 
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.size = 0.5,
               fill = 'white',
               color = 'black',
               size = 0.5
               )+
  theme_bw()+ 
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.y=element_text(size = 15), 
        axis.text.y=element_text(size=15, color = 'black'),
        axis.text.x=element_text(size=14, color = 'black', angle = 45, vjust = 0.6),
        legend.title=element_text(size = 15),
        legend.text=element_text(size = 12), 
        panel.border = element_rect(size = 2, color = 'black'), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  labs(title = '', x = '', y = 'Exchange flux (fmol / cell / h)')+
  guides(fill=guide_legend(title = NULL))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.3)

P1 <- gg.gap(plot = P,
            segments=c(40, 50),
            tick_width = c(20, 40),
            rel_heights = c(0.25, 0.0, 0.08),
            ylim=c(-40, 100))
P1

ggsave('CORE of metabolites/CORE boxplot1.jpg', P1, width = 7, height = 4,
       limitsize = FALSE, dpi = 1000)

### Gln
CORE_only_gln <- CORE_data_long %>%
  filter(Metabolite == 'Gln')

P2 <- ggplot(CORE_only_gln, 
            aes(x=Metabolite, y=Flux, fill=Metabolite)) + 
  # geom_violin(width=1, trim=TRUE,fill="#adb5bd",color="white") + 
  geom_boxplot(width=0.3,
               position=position_dodge(0.9),
               outlier.size = 0.5,
               fill = 'white',
               color = 'black',
               size = 0.5
  )+
  theme_bw()+ 
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.y=element_text(size = 11), 
        axis.text.y=element_text(size=11, color = 'black'),
        axis.text.x=element_text(size=10, color = 'black', angle = 0, vjust = 0.6),
        legend.title=element_text(size = 15),
        legend.text=element_text(size = 12), 
        panel.border = element_blank(),
        axis.line.x.bottom = element_line(color = 'black', size = 0.3),  # 只画下边 x 轴
        axis.line.y.left   = element_line(color = 'black', size = 0.3),  # 只画左边 y 轴
        axis.line.x.top    = element_blank(),   # 不画上边
        axis.line.y.right  = element_blank(),   # 不画右边
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  labs(title = '', x = '', y = 'Exchange flux (fmol / cell / h)')+
  guides(fill=guide_legend(title = NULL))

P2

ggsave('CORE of metabolites/CORE boxplot2.jpg', P2, width = 2, height = 2,
       limitsize = FALSE, dpi = 1000)






