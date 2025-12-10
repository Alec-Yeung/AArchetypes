rm(list=ls())
library(dplyr)
library(openxlsx)

# data <- read.xlsx('Combined_SubcellularLocation_in_pathways.xlsx')
# 
# location_data <- data %>%
#   dplyr::select(Cluster, Subcellular_location)
# location_data <- location_data %>% filter(!is.na(Subcellular_location))
# location_data$Subcellular_location <- tolower(location_data$Subcellular_location)
# 
# location_data$Subcellular_location <- ifelse(grepl("membrane", location_data$Subcellular_location), "Membrane protein", "Other protein")
# save(location_data, file = 'Location_data.Rdata')

### Hypergeometric Test 
load('Location_data.Rdata')

m <- sum(location_data$Subcellular_location == 'Membrane protein')
n <- sum(location_data$Subcellular_location == 'Other protein')

## A
k_A <- sum(location_data$Cluster == 'A')
q_A_membrane <- sum(location_data$Cluster == "A" & location_data$Subcellular_location == "Membrane protein")
P_A_membrane <- phyper(q_A_membrane-1, m, n, k_A, lower.tail = FALSE)

## B
k_B <- sum(location_data$Cluster == 'B')
q_B_membrane <- sum(location_data$Cluster == "B" & location_data$Subcellular_location == "Membrane protein")
P_B_membrane <- phyper(q_B_membrane-1, m, n, k_B, lower.tail = FALSE)

## C
k_C <- sum(location_data$Cluster == 'C')
q_C_membrane <- sum(location_data$Cluster == "C" & location_data$Subcellular_location == "Membrane protein")
P_C_membrane <- phyper(q_C_membrane-1, m, n, k_C, lower.tail = FALSE)

## D
k_D <- sum(location_data$Cluster == 'D')
q_D_membrane <- sum(location_data$Cluster == "D" & location_data$Subcellular_location == "Membrane protein")
P_D_membrane <- phyper(q_D_membrane-1, m, n, k_D, lower.tail = FALSE)

### Bar Plot
df_count <- location_data %>%
  group_by(Cluster, Subcellular_location) %>%
  summarise(Count = n(), .groups = 'drop')

stars <- c("A" = "ns", "B" = "ns", "C" = "****", "D" = "****")

df_count$stars <- ifelse(df_count$Subcellular_location == "Membrane protein", 
                         stars[df_count$Cluster], 
                         NA)  

P <- ggplot(df_count, aes(x = Cluster, y = Count, fill = Subcellular_location)) + 
  geom_bar(stat = 'identity', position = 'dodge', width = 0.6) + 
  geom_text(aes(label = Count), size = 3, position = position_dodge(0.6), vjust = -0.5, color = "black") +
  geom_text(aes(label = stars), size = 3, position = position_dodge(0.6), vjust = -3, color = "#cf5149") +
  labs(x = '', y = 'Count') +
  scale_fill_manual(values = c("Membrane protein" = "#6c757d", "Other protein" = "#d9d9d9")) +
  theme_bw() +
  ggtitle('Distribution of Membrane Proteins and Other Proteins by AArchetypes') +
  guides(fill = guide_legend(title = NULL)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'),
        legend.text = element_text(size = 13, color = 'black'),
        axis.title.x = element_text(size = 12, color = 'black'),
        axis.title.y = element_text(size = 13, color = 'black'),
        panel.border = element_rect(size = 1, color = 'black'), 
        plot.title = element_text(size = 14, hjust = 0.2)
  )
P

ggsave('Distribution of Membrane Proteins and Other Proteins by AArchetypes.jpg',P,width = 6, height = 5, limitsize = FALSE, dpi = 1000)


