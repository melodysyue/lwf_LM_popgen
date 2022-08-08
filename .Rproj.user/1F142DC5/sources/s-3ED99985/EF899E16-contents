library(tidyverse)

rm(list=ls())


snps <- read.table("./data/plink/lwf.pop.bim", header = F)

allpops <- read.table("./data/neutral_outlier/lwf.pcadapt.ol.qval0.01.allpops.snps", header=T)
west <- read.table("./data/neutral_outlier/lwf.pcadapt.ol.qval0.01.west.snps", header=T)
east <- read.table("./data/neutral_outlier/lwf.pcadapt.ol.qval0.01.east.snps", header=T)
ne <- read.table("./data/neutral_outlier/lwf.pcadapt.ol.qval0.01.ne.snps", header=T)

ol <- rbind(allpops, west,east,ne) %>% 
  unique() %>% 
  pull(x)

length(ol)

write.table(ol, "./data/neutral_outlier/lwf.ol.qval0.01.snps", quote=F, row.names = F)
