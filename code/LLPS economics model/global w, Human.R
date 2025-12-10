rm(list=ls())
library(openxlsx)
library(dplyr)

###
proteome_LLPS_2_AA <- read.xlsx('Proteome_LLPS_2_AA, Human.xlsx')

cols <- proteome_LLPS_2_AA[, 3:22]
col_means <- colMeans(cols, na.rm = TRUE)
col_means
sum(col_means)

###
LLPS_AA_cost <- c(K=0,A=0.30,Q=11.56,P=18.46,G=22.50,E=24.45,S=110.82)
other_AA_cost <- c(H=0,I=0,L=0,M=0,F=0,T=0,W=0,V=0,D=1.14,
                      Y=3.99,C=8.66,R=25.99,N=28.80)

common_AA <- intersect(names(col_means), names(LLPS_AA_cost))
other_AA <- setdiff(names(col_means), common_AA)

LLPS_cost <- sum(col_means[common_AA] * LLPS_AA_cost[common_AA]) / sum(col_means[common_AA])
other_cost <- sum(col_means[other_AA] * other_AA_cost[other_AA]) / sum(col_means[other_AA])

w <- LLPS_cost / other_cost






