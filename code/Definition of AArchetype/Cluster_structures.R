rm(list=ls())
library(DECIPHER)
library(openxlsx)

data <- read.table("3Di_sequence.txt", header = FALSE)
colnames(data) <- "3Di_sequence"

Alphafold_pdb_files <- list.files("Alphafold_pdb_files",pattern = "\\.pdb$",
                                  full.names = FALSE)
protein_ID <- sub("\\.pdb$", "", Alphafold_pdb_files)
df_protein_ID <- data.frame(Protein_ID = protein_ID)

rownames(data) <- df_protein_ID$Protein_ID

sequences <- data$'3Di_sequence'
sequences_set <- AAStringSet(sequences)
names(sequences_set) <- df_protein_ID$Protein_ID

set.seed(123)

clusters1 <- Clusterize(sequences_set, cutoff=seq(0.8, 0, -0.2),
                        singleLinkage = TRUE)
apply(clusters1, 2, max) # number of clusters per cutoff
apply(clusters1, 2, function(x) which.max(table(x))) # max sizes

# write.xlsx(clusters1, 'Cluster_result.xlsx', rowNames = TRUE)


aligned_seqs1 <- AlignSeqs(sequences_set, verbose=TRUE)
d1 <- DistanceMatrix(aligned_seqs1, verbose=TRUE)

# save(d1, file = '3Di_DistanceMatrix.Rdata')


