library(tidyverse)

rm(list=ls())

snps <- read.table("./data/plink/lwf.pop.bim", header=F)
snps <- snps %>% 
  select(chr=V1, pos=V4, snp=V2)

allpops <- read.table("./data/neutral_outlier/lwf.pcadapt.ol.qval0.01.allpops.snps", header=T)
west <- read.table("./data/neutral_outlier/lwf.pcadapt.ol.qval0.01.west.snps", header=T)
east <- read.table("./data/neutral_outlier/lwf.pcadapt.ol.qval0.01.east.snps", header=T)
ne <- read.table("./data/neutral_outlier/lwf.pcadapt.ol.qval0.01.ne.snps", header=T)

ol <- rbind(allpops, west,east,ne) %>% 
  unique() %>% 
  pull(x) %>% 
  as.character()



write.table(ol, "./data/neutral_outlier/lwf.ol.qval0.01.snps", quote=F, row.names = F)



#regions
r.all <- read.table("./data/islands/sigRegions_allpops.txt", header=T)
r.west <- read.table("./data/islands/sigRegions_west.txt", header=T)
r.east <- read.table("./data/islands/sigRegions_east.txt", header=T)
r.ne <- read.table("./data/islands/sigRegions_ne.txt", header=T)

r.all$data <- "All populations"
r.west$data <- "Northwest"
r.east$data <- "East"
r.ne$data <- "Northeast"


dd.r <- rbind(r.all, r.west, r.east, r.ne)

regions <- dd.r %>% 
  group_by(chromosome) %>% 
  summarize(start=min(start), end=max(end))

#region ol SNPs

regions.ol.snps <- NULL

for (i in 1:nrow(regions)){
  c=regions[i,]$chromosome
  print(paste0("Processing chr", c))
  start=regions[i,]$start
  end=regions[i,]$end
  
  temp <- snps %>% 
    filter(snp %in% ol) %>% 
    filter(chr == c) %>% 
    filter(pos>=start & pos<=end)
  
  temp.snps <- temp %>% 
    pull(snp)
  
  print(length(temp.snps))
  
  write.table(temp.snps, paste0("./data/islands/chr", c, "_region_ol_snps.txt"), quote=F, row.names = F)
  
  
regions.ol.snps <- rbind(temp, regions.ol.snps)
  
}

table(regions.ol.snps$chr)

write.table(regions.ol.snps, "./data/islands/sigRegions_ol_snps.txt", quote=F, row.names = F)


head(regions.ol.snps)
regions.size <- regions.ol.snps %>% 
  group_by(chr) %>% 
  summarize(start=min(pos), end=max(pos))

write.table(regions.size, "./data/islands/sigRegions_start_end.txt", quote=F, row.names = F)

