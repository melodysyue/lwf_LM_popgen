library(tidyverse)
library(gghighlight)
library(evobiR)

rm(list=ls())

#loci id and fst

#comparison
dd <- read.table("./data/neutral_outlier/lwf.pcadapt.east.qvalues.txt", header=T)
chr <- seq(1,39,1)


#sliding window with a window size of 250 SNPs and a step of 50 SNPs;
#significant windows are those with >= 10 outlier SNPs (qval<0.01) 

window=250
step=50
cutoff=10


result <- NULL
qval.cutoff <- 0.01

for (c in chr){
  print(paste0("Processing chr ", c))
    dd_chr <-dd %>%
      filter(chr == c)
    
    total <- nrow(dd_chr)
    spots <- seq(from = 1, to = (total - window + 1), by = step)
    
    for(i in 1:length(spots)){
      
      start <- dd_chr[spots[i], "pos"]
      end <- dd_chr[spots[i]+window-1, "pos"]
      w.ol <- dd_chr[spots[i]:(spots[i] + window - 1),] %>% 
        filter(qval<qval.cutoff) 
      
      temp <- c(index=i, chromosome=c, start=start, end=end, ol=nrow(w.ol))
      
      result <- rbind(result, temp) %>% 
        as.data.frame()
    }
    
  }


result.f <- result %>% 
    filter(ol >= cutoff) 


##combine consecutive windows into regions
result.combine <- result.f %>% 
  select(index, chromosome, start, end) %>% 
  unique() %>% 
  arrange(chromosome, index)

y=result.combine$index
g <- cumsum(c(1, abs(y[-length(y)] - y[-1]) > 1)) #whether each neighboring pair is dissimilar or not with 1 as initial value.

result.combine$region <- g


result.combine %>% 
  group_by(region) %>% 
  mutate(region_start=min(start), region_end=max(end)) %>% 
  mutate(region_size=region_end-region_start) %>% 
  select(chromosome, region, region_start, region_end, region_size) %>% 
  unique()

write.table(result.combine, "./data/islands/sigRegions_ne.txt", quote=F, row.names = F)
