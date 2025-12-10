rm(list=ls())
library(data.table)
library(openxlsx)
library(dplyr)

###
### Hypergeometric Test 
load('Combined_uniprotID+LLPS_in_pathways.Rdata')

data_AB <- merge_df2 %>%
  filter(Cluster %in% c('A', 'B'))
data_CD <- merge_df2 %>%
  filter(Cluster %in% c('C', 'D'))

m_AB <- sum(data_AB$Is_LLPS_protein == 1)
n_AB <- sum(data_AB$Is_LLPS_protein == 0)

m_CD <- sum(data_CD$Is_LLPS_protein == 1)
n_CD <- sum(data_CD$Is_LLPS_protein == 0)

## AB-A
k_A <- sum(data_AB$Cluster == 'A')
q_A <- sum(data_AB$Cluster == "A" & data_AB$Is_LLPS_protein == 1)
P_A <- phyper(q_A-1, m_AB, n_AB, k_A, lower.tail = FALSE)

## AB-B
k_B <- sum(data_AB$Cluster == 'B')
q_B <- sum(data_AB$Cluster == "B" & data_AB$Is_LLPS_protein == 1)
P_B <- phyper(q_B-1, m_AB, n_AB, k_B, lower.tail = FALSE)

## CD-C
k_C <- sum(data_CD$Cluster == 'C')
q_C <- sum(data_CD$Cluster == "C" & data_CD$Is_LLPS_protein == 1)
P_C <- phyper(q_C-1, m_CD, n_CD, k_C, lower.tail = FALSE)

## CD-D
k_D <- sum(data_CD$Cluster == 'D')
q_D <- sum(data_CD$Cluster == "D" & data_CD$Is_LLPS_protein == 1)
P_D <- phyper(q_D-1, m_CD, n_CD, k_D, lower.tail = FALSE)

### AB
df_count_AB <- data_AB %>%
  group_by(Cluster, Is_LLPS_protein) %>%
  summarise(Count = n(), .groups = 'drop')

df_count_AB$Is_LLPS_protein <- factor(df_count_AB$Is_LLPS_protein,
                                   levels = c(1, 0),
                                   labels = c("LLPS Protein", "Not-LLPS Protein"))

P <- ggplot(df_count_AB, aes(x = Cluster, y = Count, fill = Is_LLPS_protein)) + 
  geom_bar(stat = 'identity', position = 'dodge', width = 0.6) + 
  geom_text(aes(label = Count), size = 3, position = position_dodge(0.6), vjust = -0.5, color = "black") +
  # geom_text(aes(label = stars), size = 3, position = position_dodge(0.6), vjust = -3, color = "#cf5149") +  
  labs(x = '', y = 'Count') +
  scale_fill_manual(values = c("LLPS Protein" = "#6c757d", "Not-LLPS Protein" = "#d9d9d9")) +
  theme_bw() +
  ggtitle('') +
  guides(fill = guide_legend(title = NULL)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 15, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'),
        legend.text = element_text(size = 11, color = 'black'),
        axis.title.x = element_text(size = 12, color = 'black'),
        axis.title.y = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 12, color = 'black'),
        plot.title = element_text(size = 15, hjust = 0.3),
        panel.border = element_rect(size = 1, color = 'black')
  )
P

### CD
df_count_CD <- data_CD %>%
  group_by(Cluster, Is_LLPS_protein) %>%
  summarise(Count = n(), .groups = 'drop')

df_count_CD$Is_LLPS_protein <- factor(df_count_CD$Is_LLPS_protein,
                                   levels = c(1, 0),
                                   labels = c("LLPS Protein", "Not-LLPS Protein"))

stars_CD <- c("C" = "ns", "D" = "****")

df_count_CD$stars <- ifelse(
  df_count_CD$Is_LLPS_protein == 'Not-LLPS Protein', NA,  # 如果 Is_LLPS_protein 为 Not-LLPS Protein，stars 为 NA
  ifelse(
    df_count_CD$Cluster == "D", "****",  # 如果 Is_LLPS_protein 为 1 且 Cluster 为 D，stars 为 "***"
    "ns"                       # 如果 Is_LLPS_protein 为 1 且 Cluster 不为 D，stars 为 "ns"
  )
)

P <- ggplot(df_count_CD, aes(x = Cluster, y = Count, fill = Is_LLPS_protein)) + 
  geom_bar(stat = 'identity', position = 'dodge', width = 0.6) + 
  geom_text(aes(label = Count), size = 3, position = position_dodge(0.6), vjust = -0.5, color = "black") +
  geom_text(aes(label = stars), size = 3, position = position_dodge(0.6), vjust = -3, color = "#cf5149") +  
  labs(x = '', y = 'Count') +
  scale_fill_manual(values = c("LLPS Protein" = "#6c757d", "Not-LLPS Protein" = "#d9d9d9")) +
  theme_bw() +
  ggtitle('') +
  guides(fill = guide_legend(title = NULL)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 15, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'),
        legend.text = element_text(size = 11, color = 'black'),
        axis.title.x = element_text(size = 12, color = 'black'),
        axis.title.y = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 12, color = 'black'),
        plot.title = element_text(size = 15, hjust = 0.3),
        panel.border = element_rect(size = 1, color = 'black')
  )
P














