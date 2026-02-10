rm(list=ls())
library(openxlsx)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(fgsea)
library(dplyr)
library(ggplot2)
library(data.table)

###
df_cost <- read.xlsx('Human_proteome_cost.xlsx')
df_cost$Protein_cost <- as.numeric(df_cost$Protein_cost)
df_cost$Energy_cost <- as.numeric(df_cost$Energy_cost)

df_cost_decreasing <- df_cost[order(df_cost$Protein_cost, decreasing = TRUE), ]

###
df_cost_decreasing$EntrezID <- mapIds(org.Hs.eg.db,
                                      keys = df_cost_decreasing$UniprotID,
                                      column = "ENTREZID",     
                                      keytype = "UNIPROT",     
                                      multiVals = "first")     

###
df_cost_decreasing <- df_cost_decreasing %>% filter(!is.na(EntrezID))

###
genelist <- df_cost_decreasing$Protein_cost
names(genelist) <- as.character(df_cost_decreasing$EntrezID)

genelist <- genelist[!duplicated(names(genelist))]

head(genelist)

###
IDR_data <- read.csv("D:\\Desktop\\Sequence_complexity\\IDR_data.csv")

IDR_data1 <- IDR_data %>% dplyr::select(UniProtID)
colnames(IDR_data1)[1] <- 'UniprotID'

IDR_data1$UniprotID <- sub("-\\d+$", "", IDR_data1$UniprotID) 
IDR_data1$UniprotID <- trimws(IDR_data1$UniprotID)             

IDR_data2 <- IDR_data1
IDR_data2$EntrezID <- mapIds(org.Hs.eg.db,
                             keys = IDR_data2$UniprotID,
                             column = "ENTREZID",     
                             keytype = "UNIPROT",     
                             multiVals = "first")     
IDR_data2 <- unique(IDR_data2)
IDR_data2 <- IDR_data2 %>% filter(!is.na(EntrezID))

IDR_genes <- IDR_data2$EntrezID
IDR_genes <- unique(IDR_genes)

###
LLPS_data <- fread("D:\\Desktop\\LLPS\\MobiDB.txt")

LLPS_human <- LLPS_data %>%
  filter(grepl('Hos', `LLPS ID`)) %>%
  dplyr::select(`UniProt ID`)

LLPS_human$EntrezID <- mapIds(org.Hs.eg.db,
                              keys = LLPS_human$`UniProt ID`,
                              column = "ENTREZID",   
                              keytype = "UNIPROT",    
                              multiVals = "first")     
LLPS_human <- LLPS_human %>% filter(!is.na(EntrezID))

LLPS_genes <- LLPS_human$EntrezID

###
kegg <- gmtPathways("c2.cp.v2025.1.Hs.entrez.gmt")

pathway_all <- c(kegg, list('IDR_pathway' = IDR_genes,
                            'LLPS_pathway' = LLPS_genes))

###
res <- fgsea(pathways = pathway_all,
             stats = genelist,
             minSize = 10, maxSize = 20000, nperm = 10000)

###
IDR_NES <- 6.6
LLPS_NES <- 3.8

vlines <- data.frame(
  x = c(IDR_NES, LLPS_NES), 
  y = 0,
  yend = 100)

P <- ggplot(res, aes(x = NES))+
  geom_histogram(binwidth = 0.1, fill = "#b7b9b8", color = "#b7b9b8")+
  labs(x = 'NES',
       y = 'Number of KEGG pathways',
       title = 'Distribution of NES')+
  theme_minimal()+
  geom_segment(data = vlines, 
               aes(x = x, xend = x, y = y, yend = yend), inherit.aes = FALSE,
               color = "#cf5149", size = 1, linetype = "solid")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(size = 17, hjust = 0.5),
        axis.title.y = element_text(size = 16), 
        axis.title.x = element_text(size = 16),
        axis.text.y=element_text(size=13, color = 'black'),
        axis.text.x=element_text(size=13, color = 'black'))
P

ggsave('Distribution of NES - proteome.jpg',P,
       width = 5, height = 5, limitsize = FALSE, dpi = 1000)

###
P <- plotEnrichment(pathway_all[['IDR_pathway']], genelist) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(size = 17, hjust = 0.5),
        axis.title.y = element_text(size = 16), 
        axis.title.x = element_text(size = 16),
        axis.text.y=element_text(size=12, color = 'black'),
        axis.text.x=element_text(size=12, color = 'black'))

P

###
P <- plotEnrichment(pathway_all[['LLPS_pathway']], genelist) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(size = 17, hjust = 0.5),
        axis.title.y = element_text(size = 16), 
        axis.title.x = element_text(size = 16),
        axis.text.y=element_text(size=13, color = 'black'),
        axis.text.x=element_text(size=13, color = 'black'))
P





#######################################
###
ekegg <- gseKEGG(
  geneList     = genelist,
  organism     = "hsa",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 1,
  verbose      = FALSE
)

ekegg_res <- as.data.frame(ekegg) |> dplyr::select(ID, Description, NES, pvalue, p.adjust, setSize)

ridgeplot(ekegg, showCategory = 20)

###
df_pathway <- ekegg_res %>%
  mutate(group = ifelse(NES > 0, "+", "-")) %>%  # 以A列正负分组
  group_by(group) %>%
  slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>% # 取B列最小的10行
  ungroup() %>%
  dplyr::select(NES, p.adjust, Description)   # 保留ABC三列

###
# 确保 Description 按 NES 顺序排列（方便从上到下按数值显示）
df_pathway$Description <- factor(df_pathway$Description, levels = df_pathway$Description[order(df_pathway$NES)])

###
df_pathway$`-log(p.adjust)` <- -log10(df_pathway$p.adjust)

###
P <- ggplot(df_pathway, aes(x = NES, y = Description, fill = `-log(p.adjust)`)) +
  geom_col() +
  scale_fill_gradient(low = "#9ecafe", high = "#fecee6",
                      name = "-log(p.adjust)") +
  geom_vline(xintercept = 0, color = "white", linetype = "solid") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 13, color = 'black'),
    axis.text.x = element_text(size = 12, color = 'black'),
    legend.position = "right"
  ) +
  labs(
    x = "NES",
    title = "GSEA NES by KEGG pathway"
  )
P

ggsave('GSEA NES by pathway-proteome.jpg',P, width = 10, height = 7, limitsize = FALSE, dpi = 1000)











