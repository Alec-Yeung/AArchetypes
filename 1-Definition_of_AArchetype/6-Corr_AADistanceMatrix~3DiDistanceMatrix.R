rm(list = ls())
library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggrastr)

### 
load("D:\\Desktop\\LLPS\\Combined_uniprotID+LLPS_in_pathways.Rdata")
Combined_genes_in_pathways <- read.xlsx("D:\\Desktop\\AA_composition\\Combined_genes_in_pathways.xlsx")
load('3Di_DistanceMatrix.Rdata')

uniprotid <- merge_df2 %>%
  dplyr::select(Gene, uniprotswissprot)
uniprotid <- unique(uniprotid)

AA_compostion <- Combined_genes_in_pathways %>%
  dplyr::select(Gene, A:Y)
AA_compostion <- unique(AA_compostion)
AA_compostion$Gene <- as.character(AA_compostion$Gene)

proteinID <- rownames(d1)

### 
merge_df3 <- uniprotid %>%
  left_join(AA_compostion, by = 'Gene')

###
merge_df3 <- merge_df3[merge_df3$uniprotswissprot %in% proteinID, ]

merge_df3 <- merge_df3[match(proteinID, merge_df3$uniprotswissprot), ]

merge_df3 <- na.omit(merge_df3)

###
AA_compostion1 <- merge_df3
AA_compostion1$Gene <- NULL
rownames(AA_compostion1) <- AA_compostion1[, 1]   
AA_compostion1 <- AA_compostion1[, -1]

d2 <- dist(AA_compostion1, method = "euclidean")
d2 <- as.matrix(d2)

###
common_proteins <- rownames(d2)

d3 <- d1[common_proteins, common_proteins]

all(rownames(d2) == rownames(d3)) # TRUE
all(colnames(d2) == colnames(d3)) # TRUE

###
d2_vector <- d2[upper.tri(d2)]
d3_vector <- d3[upper.tri(d3)]

sum(is.na(d2_vector)) # 0
sum(is.na(d3_vector)) # 1985

cor_test_result <- cor.test(d2_vector, d3_vector, method = "pearson")
print(cor_test_result)

###
data <- data.frame(AACompositionDistance = d2_vector, `3DiSequenceDistance` = d3_vector)

correlation <- cor_test_result$estimate 

p <- cor_test_result$p.value 

###
P <- ggplot(data, aes(x=`X3DiSequenceDistance`, y=AACompositionDistance)) + 
  # geom_point(color = '#b7b9b8', na.rm = TRUE)+ 
  geom_point_rast(color = '#71b0fd', size = 0.3, na.rm = TRUE) + # b7b9b8
  theme_minimal()+
  labs(title = '',
       x = 'Structural Distance', 
       y = 'AA Composition Distance')+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(size = 17, hjust = 0.5),
        axis.title.y = element_text(size = 16), 
        axis.title.x = element_text(size = 16),
        axis.text.y=element_text(size=13, color = 'black'),
        axis.text.x=element_text(size=13, color = 'black'))+
  # annotate("text", x = Inf, y = Inf, label = paste("Pearson's r =", round(correlation, 3)),
  #          hjust = 2, vjust = 2, size = 4, color = "black")+
  stat_density_2d(aes(fill = after_stat(level)), 
                  geom = "polygon", # raster
                  contour_var = "density",
                  bins = 10,
                  colour = NA,
                  alpha = 1) +
  scale_fill_gradientn(colors = c("#71b0fd", "#9ecafe", "#c4dffc", # "#a8d8ca", "#bcdd9d", "#e9e944"
                                  "#ffeef8", "#fecee6", "#fdc0e1" # "#fbd723", "#fa9f30", "#f37331"
                                  ), # "#ed3536"
                       name = "")+
  guides(fill = guide_colorbar(title = "Density",
                               title.position = "top",
                               barwidth = 1,
                               barheight = 5))
P

ggsave('Correlation between AACompositionDistance and 3DiSequenceDistance.tiff',
       P, width = 5, height = 5, limitsize = FALSE, dpi = 1000, compression = "lzw")







