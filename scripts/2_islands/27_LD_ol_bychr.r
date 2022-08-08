library(tidyverse)
library(patchwork)

rm(list=ls())

###LD heatmap###
ld.lists <- list.files(path = "./data/ld_r2/", pattern = "*.ld", full.names=TRUE) %>% 
  str_sort(numeric = TRUE)

ld_df <- ld.lists %>% 
  map_df(read.table, header=T) 

chrs <- list.files(path = "./data/ld_r2/", pattern = "*.ld", full.names=TRUE) %>% 
  str_sort(numeric=TRUE) %>% 
  basename() %>% 
  str_remove_all(., '.ld') %>% 
  as.numeric()

ld.lists
chrs
n <- length(ld.lists)
table(ld_df$CHR_A)

ol <- read.table("./data/islands/sigRegions_ol_snps.txt", header=T)


ld_df_ol <- NULL

for(i in 1:n){
  c=chrs[i]
  print(paste0("Processing chromosome ", c))
  ld_chr <- ld_df %>% 
    filter(CHR_A==c)
  
  ol_chr <- ol %>% 
    filter(chr==c)
  
  ld_chr <- ld_chr %>% 
    mutate(pairs = ifelse(SNP_A %in% ol_chr$snp & SNP_B %in% ol_chr$snp, "Outliers", "Background"))
  ld_chr$pairs <- as.factor(ld_chr$pairs)
  ld_chr$pairs <- factor(ld_chr$pairs, levels=c("Outliers", "Background"))
  
  ld_df_ol <- rbind(ld_df_ol, ld_chr)

}


ld_df_ol %>% 
  group_by(CHR_A, pairs) %>%
  summarize(n=n()) %>% 
  spread(pairs,n)
  

ld_df_ol$CHR_A <- as.factor(ld_df_ol$CHR_A)
levels(ld_df_ol$CHR_A)
ld_df_ol$CHR_A <- factor(ld_df_ol$CHR_A, 
                  levels=c("4", "7", "10", "11", "18", "20"),
                  labels=c("Chromosome 4", "Chromosome 7", "Chromosome 10", "Chromosome 11", "Chromosome 18", "Chromosome 20"))


###plot
###outliers vs background
p_ld <- ld_df_ol %>% 
  ggplot(aes(x=pairs, y=R2, fill=pairs))+
  geom_boxplot(alpha=0.5)+
  ylab(expression(italic(r^2)))+
  scale_fill_manual(values=c("red", "dimgray"), labels=c("Pcadapt outliers", "Background"))+
  #scale_fill_manual(values=c("red", "dimgray"), labels=c(expression(italic(pcadapt)~' outliers'), "Background"))+
  theme_classic(base_size=20)+
  theme(axis.title.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=1),
        legend.position = "bottom",
        legend.title = element_blank())+
  facet_wrap(~CHR_A)


pdf("./figures/figS6_LD_ol_background_compare.pdf", width = 12, height = 9)
p_ld
dev.off()



