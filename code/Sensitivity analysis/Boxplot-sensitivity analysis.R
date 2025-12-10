rm(list = ls())
library(openxlsx)
library(ggplot2)
library(tidyr)
library(gg.gap)

###
data <- read.xlsx('CostSensitivity_of_pathways_to_AAs_e.xlsx') ##

plot_data <- data %>%
  pivot_longer(
    cols = -Pathway, 
    names_to = "AA", 
    values_to = "SensitivityScore" 
  )

plot_data$Pathway <- NULL

###
median_score <- plot_data %>%
  group_by(AA) %>%
  summarise(median_score = median(SensitivityScore, na.rm = TRUE))

ordered_levels <- plot_data %>%
  group_by(AA) %>%
  summarise(median_score = median(SensitivityScore, na.rm = TRUE)) %>%
  arrange(median_score) %>% 
  pull(AA) 

plot_data <- plot_data %>%
  mutate(AA = factor(AA, levels = ordered_levels))

###
P <- ggplot(plot_data, 
            aes(x=AA, y=SensitivityScore, fill=AA)) + 
  geom_boxplot(size=0.1,width=0.5,position=position_dodge(0.9),outlier.colour = NA, fill = 'white')+
  geom_jitter(size = 0.5, width = 0.1, alpha = 0.2, color = "#cf5149") +  
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.3) + 
  stat_summary(fun = median, 
               geom = "text", 
               aes(label = round(..y.., 2)), 
               angle = 30,
               vjust = -1.2,
               hjust = 0.1,
               size = 4, 
               color = "black") +
  theme_bw()+ 
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.y=element_text(size = 15), 
        axis.text.y=element_text(size=11, color = 'black'),
        axis.text.x=element_text(size=10.5, color = 'black', angle = 40, hjust = 0.7),
        legend.title=element_text(size = 15),
        legend.text=element_text(size = 12), 
        legend.position = "none",
        panel.border = element_rect(size = 1, color = 'black'), 
        # axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  labs(title = 'Sensitivity Scores of AAs', x = '', y = '')+
  guides(fill=guide_legend(title = NULL))
P


P1 <- gg.gap(plot = P,
             segments=c(5, 12),
             tick_width = c(2, 2),
             rel_heights = c(0.03, 0.0, 0.01),
             ylim=c(-2, 15)) 

P1

ggsave('SensitivityScores of AAs2.png', P, width = 5, height = 4,
       limitsize = FALSE, dpi = 1000)






