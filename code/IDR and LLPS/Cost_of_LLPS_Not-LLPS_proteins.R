rm(list = ls())
library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggpubr)

### LLPS data
load('Combined_uniprotID+LLPS_in_pathways.Rdata')

LLPS_data <- merge_df2 %>%
  dplyr::select(Is_LLPS_protein, Gene)

LLPS_data <- unique(LLPS_data)

### Cost data
cost_data <- read.xlsx("D:\\Desktop\\Cost\\Combined_costs_of_proteins.xlsx")

cost_data <- cost_data %>%
  dplyr::select(Gene, Protein_cost, Energy_cost)

cost_data <- unique(cost_data)
cost_data$Gene <- as.character(cost_data$Gene)

### merge LLPS_data & cost_data
merge_df3 <- LLPS_data %>%
  left_join(cost_data, by = 'Gene')

### Violin Plot: Protein cost
violin_data1 <- merge_df3 %>%
  dplyr::select(Is_LLPS_protein, Protein_cost)

violin_data1 <- violin_data1 %>%
  mutate(Is_LLPS_protein = ifelse(Is_LLPS_protein == 1, "LLPS Protein", ifelse(Is_LLPS_protein == 0, "Not-LLPS Protein", Is_LLPS_protein)))
violin_data1$Protein_cost <- as.numeric(violin_data1$Protein_cost)

P <- ggplot(violin_data1, 
            aes(x=Is_LLPS_protein, y=Protein_cost, fill=Is_LLPS_protein)) + 
  geom_violin(width=0.8, trim=TRUE,fill="#adb5bd",color="white", drop = FALSE) + 
  geom_boxplot(width=0.1,position=position_dodge(0.9),outlier.colour = NA, fill = 'white')+
  theme_bw()+ 
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.y=element_text(size = 17), 
        axis.text.y=element_text(size=15, color = 'black'),
        axis.text.x=element_text(size=16, color = 'black'),
        legend.title=element_text(size = 15),
        legend.text=element_text(size = 12), 
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  labs(title = 'Protein Cost of LLPS/Not-LLPS Protein', x = '', y = '')+
  guides(fill=guide_legend(title = NULL))+
  stat_compare_means(comparisons = list(c("LLPS Protein", "Not-LLPS Protein")),
                     method = "wilcox.test", label = "p.signif", size = 8, bracket.size = 0.5, vjust = 0.5,
                     tip.length = 0.02,
                     label.y = c(22)
  )
P     

ggsave('Protein Cost of LLPS and Not-LLPS Protein.jpg', P, width = 5, height = 5, limitsize = FALSE, dpi = 1000)

### Violin Plot: Energy cost
violin_data2 <- merge_df3 %>%
  dplyr::select(Is_LLPS_protein, Energy_cost)

violin_data2 <- violin_data2 %>%
  mutate(Is_LLPS_protein = ifelse(Is_LLPS_protein == 1, "LLPS Protein", ifelse(Is_LLPS_protein == 0, "Not-LLPS Protein", Is_LLPS_protein)))
violin_data2$Energy_cost <- as.numeric(violin_data2$Energy_cost)

P <- ggplot(violin_data2, 
            aes(x=Is_LLPS_protein, y=Energy_cost, fill=Is_LLPS_protein)) + 
  geom_violin(width=0.8, trim=TRUE,fill="#adb5bd",color="white", drop = FALSE) + 
  geom_boxplot(width=0.1,position=position_dodge(0.9),outlier.colour = NA, fill = 'white')+
  theme_bw()+ 
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.y=element_text(size = 17), 
        axis.text.y=element_text(size=15, color = 'black'),
        axis.text.x=element_text(size=16, color = 'black'),
        legend.title=element_text(size = 15),
        legend.text=element_text(size = 12), 
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  labs(title = 'Energy Cost of LLPS/Not-LLPS Protein', x = '', y = '')+
  guides(fill=guide_legend(title = NULL))+
  stat_compare_means(comparisons = list(c("LLPS Protein", "Not-LLPS Protein")),
                     method = "wilcox.test", label = "p.signif", size = 8, bracket.size = 0.5, vjust = 0.5,
                     tip.length = 0.02,
                     label.y = c(2.7)
  )
P     

ggsave('Energy Cost of LLPS and Not-LLPS Protein.jpg', P, width = 5, height = 5, limitsize = FALSE, dpi = 1000)





