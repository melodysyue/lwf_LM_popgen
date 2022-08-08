library(tidyverse)
library(LDheatmap)
library(stringr)

rm(list=ls())
###LD heatmap###
ld.lists <- list.files(path = "./data/ld_r2/", pattern = "*.ld", full.names=TRUE) %>% 
  str_sort(numeric = TRUE)

ld_df <- ld.lists %>% 
  lapply(read.table, header=TRUE)

chrs <- list.files(path = "./data/ld_r2/", pattern = "*.ld", full.names=TRUE) %>% 
  str_sort(numeric=TRUE) %>% 
  basename() %>% 
  str_remove_all(., '.ld')

ld.lists
chrs

n <- length(ld.lists)

region <- read.table("./data/islands/sigRegions_start_end.txt", header=T)
head(region)

#LDheatmap

for(i in 1:n){
  print(paste0("Processing chromosome ", chrs[i]))
  #input: square matrix with values between 0 - 1 inclusive.
  
  ld <- ld_df[[i]]
  ld <- ld %>%
    select(BP_A, BP_B, R2)
  
  nameVals <- sort(unique(unlist(ld[1:2])))
  nameVals <- as.character(nameVals)
  
  # construct 0 matrix of correct dimensions with row and column names
  myMat <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))

  ld$BP_A <- as.character(ld$BP_A)
  ld$BP_B <- as.character(ld$BP_B)
  # fill in the matrix with matrix indexing on row and column names
  myMat[as.matrix(ld[c("BP_A", "BP_B")])] <- ld[["R2"]]
  
  title <- paste0("Chromosome ", chrs[[i]])
  mycols <- colorRampPalette(c("lightgray", "blue", "red"))
  
  
  start <- region[i,]$start
  end <- region[i,]$end
  pos.vector <- as.numeric(nameVals)
  start.i <- min(which(pos.vector>=start))
  end.j <- max(which(pos.vector<=end))
  

  png(paste0("./figures/figs7_chr", chrs[[i]], "_ld_heatmap_highlight.png"), width=4, height=4, units="in", res=300)
  p_ld <- LDheatmap(myMat, title=title, LDmeasure = "r", add.map=FALSE, add.key=FALSE, color=rev(mycols(30)))
  LDheatmap.highlight(p_ld, i=start.i, j=end.j, col="purple")
  dev.off()
}


#LD within candidate regions

for(i in 1:n){
  print(paste0("Processing chromosome ", chrs[i]))
  #input: square matrix with values between 0 - 1 inclusive.
  
  ld <- ld_df[[i]]
  ld <- ld %>%
    select(BP_A, BP_B, R2)

  
  start <- region[i,]$start
  end <- region[i,]$end
  
  ld %>% 
    filter(BP_A >=start & BP_A <= end & BP_B >= start & BP_B <= end) %>% 
    summary() %>% 
    print()

}





