rm(list=ls())
library(ggplot2)
library(openxlsx)

###
data <- read.xlsx('Yeast_proteomic_AA_abundance+cost.xlsx')

###
result_p <- cor.test(data$Average_proteomic_compostion, data$Protein_cost,
                   method = 'spearman')
correlation_p <- round(result_p$estimate, 2)
p_value_p <- result_p$p.value

###
P <- ggplot(data, aes(x=Protein_cost, y=Average_proteomic_compostion, label = AA)) + 
  # geom_smooth(method = 'lm', color = 'black', size = 0.5) +
  geom_text(aes(label=AA), check_overlap = FALSE, hjust = 0.5, vjust = -0.5,
            size = 5, colour='black')+
  theme_minimal()+
  labs(title = '', x = 'Protein cost', 
       y = 'Proteomic abundance')+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(size = 17, hjust = 0.5),
        axis.title.y = element_text(size = 16), 
        axis.title.x = element_text(size = 16),
        axis.text.y=element_text(size=13, color = 'black'),
        axis.text.x=element_text(size=13, color = 'black'))+
  geom_point(size = 3, color = '#c55047')+ 
  annotate("text", x = Inf, y = Inf, label = paste("Spearman's ρ =", correlation_p), 
           hjust = 2, vjust = 2, size = 4, color = "black")
P




