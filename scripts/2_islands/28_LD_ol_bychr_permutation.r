library(tidyverse)
library(ggsignif) #geom_signif
library(scales)
library(patchwork)
library(ggpubr)
library(grid)


rm(list=ls())

chrs=c(4,7,10,11,18,20)

for(c in chrs) {

  print(paste0("Processing chromosome ", c))
  


ld <- read.table(paste0("./data/ld_r2/",c,".ld"), header=T)
ol <- read.table("./data/islands/sigRegions_ol_snps.txt", header=T)
ol <- ol %>% 
  filter(chr==c)


###perform permutation tests
list.files("./scripts/permutation_functions_LD/", full.names = TRUE) %>% sapply(source) %>% invisible #load all the related functions
num_permutations <- 10000

# draw 10000 samples of num.outliers random loci, take the mean, and return the ecdf and mean
null <- replicate(num_permutations, calculate.null.metrics(ld, ol))

# calculate the estimate mean null 
null.mean <- mean(null, na.rm = TRUE)

# calculate the empirical nndist for real outliers
empirical.mean <- calculate.emp.metrics(ld, ol)

#summary

perm.stat <- data.frame(
  chr = unique(ol$chr),
  num.outliers = nrow(ol),
  
  #ld
  mean.emp = empirical.mean,
  mean.null = null.mean, 
  pvalue = two_side_p(null, empirical.mean)
  )

print(perm.stat)

}



