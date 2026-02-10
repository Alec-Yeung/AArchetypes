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








########### Sequences Alignment, DistanceMatrix, Tree ###########
aligned_seqs1 <- AlignSeqs(sequences_set, verbose=TRUE)
d1 <- DistanceMatrix(aligned_seqs1, verbose=TRUE)
tree1 <- TreeLine(myDistMatrix=d1, method="UPGMA", verbose=TRUE,
                  showPlot = TRUE)

heatmap(as.matrix(clusters1), scale="column", Colv=NA, Rowv=tree1)




### R code from vignette source 'ClusteringSequences.Rnw'

###################################################
### code chunk number 1: ClusteringSequences.Rnw:51-54
###################################################
options(continue=" ")
options(width=80)
set.seed(123)


###################################################
### code chunk number 2: startup
###################################################
library(DECIPHER)


###################################################
### code chunk number 3: expr1
###################################################
# specify the path to your file of sequences:
fas <- "<<path to training FASTA file>>"
# OR use the example DNA sequences:
fas <- system.file("extdata",
                   "50S_ribosomal_protein_L2.fas",
                   package="DECIPHER")
# read the sequences into memory
dna <- readDNAStringSet(fas)
dna


###################################################
### code chunk number 4: expr2
###################################################
aa <- translate(dna)
aa
seqs <- aa # could also cluster the nucleotides
length(seqs)


###################################################
### code chunk number 5: expr3
###################################################
clusters <- Clusterize(seqs, cutoff=seq(0.5, 0, -0.1), processors=1)
class(clusters)
colnames(clusters)
str(clusters)
apply(clusters, 2, max) # number of clusters per cutoff
apply(clusters, 2, function(x) which.max(table(x))) # max sizes


###################################################
### code chunk number 6: expr4
###################################################
o <- order(clusters[[2]], width(seqs), decreasing=TRUE) # 40% cutoff
o <- o[!duplicated(clusters[[2]])]
aa[o]
dna[o]


###################################################
### code chunk number 7: expr5
###################################################
t <- table(clusters[[1]]) # select the clusters at a cutoff
t <- sort(t, decreasing=TRUE)
head(t)
w <- which(clusters[[1]] == names(t[1]))
AlignSeqs(seqs[w], verbose=FALSE)


###################################################
### code chunk number 8: expr6
###################################################
aligned_seqs <- AlignSeqs(seqs, verbose=TRUE)
d <- DistanceMatrix(aligned_seqs, verbose=TRUE)
tree <- TreeLine(myDistMatrix=d, method="UPGMA", verbose=TRUE)
heatmap(as.matrix(clusters), scale="column", Colv=NA, Rowv=tree)


###################################################
### code chunk number 9: expr7
###################################################
c1 <- Clusterize(dna, cutoff=0.2, invertCenters=TRUE, processors=1)
w <- which(c1 < 0 & !duplicated(c1))
dna[w] # select cluster representatives (negative cluster numbers)


###################################################
### code chunk number 10: expr8
###################################################
c2 <- Clusterize(dna, cutoff=0.2, singleLinkage=TRUE, processors=1)
max(abs(c1)) # center-linkage
max(c2) # single-linkage (fewer clusters, but broader clusters)


###################################################
### code chunk number 11: expr9
###################################################
genus <- sapply(strsplit(names(dna), " "), `[`, 1)
t <- table(genus, c2[[1]])
heatmap(sqrt(t), scale="none", Rowv=NA, col=hcl.colors(100))


###################################################
### code chunk number 12: expr10
###################################################
batchSize <- 2e2 # normally a large number (e.g., 1e6 or 1e7)
o <- order(width(seqs), decreasing=TRUE) # process largest to smallest
c3 <- integer(length(seqs)) # cluster numbers
repeat {
  m <- which(c3 < 0) # existing cluster representatives
  m <- m[!duplicated(c3[m])] # remove redundant sequences
  if (length(m) >= batchSize)
    stop("batchSize is too small")
  w <- head(c(m, o[c3[o] == 0L]), batchSize)
  if (!any(c3[w] == 0L)) {
    if (any(c3[-w] == 0L))
      stop("batchSize is too small")
    break # done
  }
  m <- m[match(abs(c3[-w]), abs(c3[m]))]
  c3[w] <- Clusterize(seqs[w], cutoff=0.05, invertCenters=TRUE)[[1]]
  c3[-w] <- ifelse(is.na(c3[m]), 0L, abs(c3[m]))
}
table(abs(c3)) # cluster sizes


###################################################
### code chunk number 13: expr11
###################################################
# simulate half of strands having opposite orientation
s <- sample(c(TRUE, FALSE), length(dna), replace=TRUE)
dna[s] <- reverseComplement(dna[s])

# cluster both strands at the same time
clus <- Clusterize(c(dna, reverseComplement(dna)), cutoff=0.2, processors=1)
clus <- match(clus[[1]], clus[[1]]) # renumber clusters ascending

# if needed, reorient all clustered sequences to have the same orientation
strand <- clus[seq_len(length(clus)/2)] >= clus[-seq_len(length(clus)/2)]
dna[strand] <- reverseComplement(dna[strand])

# renumber clusters across both strands and compare to original clustering
clus <- pmin(clus[seq_len(length(clus)/2)], clus[-seq_len(length(clus)/2)])
org <- match(abs(c1[[1]]), abs(c1[[1]])) # renumber original clustering
mean(clus == org) # some differences expected due to algorithm stochasticity

# verify the largest cluster is now back in the same orientation
dna[clus == which.max(tabulate(clus))]


###################################################
### code chunk number 14: expr12
###################################################
set.seed(123) # initialize the random number generator
clusters <- Clusterize(seqs, cutoff=0.1, processors=1)
set.seed(NULL) # reset the seed


###################################################
### code chunk number 15: sessinfo
###################################################
toLatex(sessionInfo(), locale=FALSE)

