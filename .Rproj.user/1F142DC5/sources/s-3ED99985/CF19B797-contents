library(tidyverse)
library(ggsci)
library(pheatmap)
library(scales)
library(RColorBrewer)

rm(list=ls())

snps <- read.table("./data/plink/lwf.pop.bim", header=F)
snps <- snps %>% 
  select(chr=V1, pos=V4, snp=V2)


list <- list.files(path="./data/afreq/", pattern="*.afreq", full.names=TRUE, recursive = TRUE) %>% 
  sort()

pops <- str_remove_all(list, "./data/afreq/plink2.") %>% 
  sort()
pops <- str_remove_all(pops, ".afreq")


out <- list %>%
  map_df(read.table, header=F) 

colnames(out) <- c("chr", "snp", "ref",'alt',"alt_freq","allele_count")

out %>% 
  select(chr, snp, ref, alt) %>% 
  unique() %>% 
  dim() #make sure the ref/ale info are consistent


regions.ol <- read.table("./data/islands/sigRegions_ol_snps.txt", header=T)

regions.ol <- regions.ol %>% 
  mutate(chr_pos=paste0(chr, "_", pos))

out <- out %>% 
  filter(snp %in% regions.ol$snp)

n=nrow(regions.ol)

order <- c( "FOXR", "RILY","LSRF", "PSTG","MENO", "WTFB","NMNL", "BBDN","MNST", 
            "NAUB", "EPFT", "CXVL", "LTTV", "EKRP", "GOOD", "IGPT", "MSKG")

af <- out %>% 
  mutate(pop=rep(pops, each=n)) %>% 
  select(-allele_count) %>% 
  spread(pop, alt_freq) %>% 
  select(chr, snp, ref, alt, order)

head(af)
head(regions.ol)
table(regions.ol$chr)

af <- left_join(af, regions.ol, by=c("snp", "chr"))
head(af)

write.csv(af, "./data/islands/sigRegions_ol_snps_afBypop.csv", quote=F, row.names = F)


af.mat <- af %>% 
  arrange(chr, pos) %>% 
  select(-chr, -snp, -pos, -ref, -alt) %>% 
  column_to_rownames("chr_pos") %>% 
  t() %>% 
  as.matrix()

##annotation

#snp annotation
chr <- regions.ol %>% 
  arrange(chr,pos) %>% 
  select(chr_pos, chr) %>% 
  remove_rownames() %>% 
  column_to_rownames("chr_pos")

chr$chr<- as.factor(chr$chr)


west <- c("FOXR", "RILY","LSRF", "PSTG","MENO", "WTFB","NMNL", "BBDN","MNST")
east <- c("NAUB", "EPFT", "CXVL", "LTTV", "EKRP", "GOOD", "IGPT", "MSKG")

pop_basin <- pops %>% 
  as.data.frame()
colnames(pop_basin) <- "pop"

pop_basin <- pop_basin %>% 
  mutate(region=case_when(
    pop%in%west ~ "Northwest",
    pop%in%east ~ "East"
  )) %>%
  column_to_rownames("pop")

pop_basin$region <- as.factor(pop_basin$region)
pop_basin$region <- factor(pop_basin$regio, levels=c("Northwest", "East"))


#specify annotation colors

region <- c("steelblue", "orange")
names(region) <- c("Northwest", "East")
anno_colors <- list(region=region)

pdf("./figures/figS5_sigRegions_ol_snps_afBypop_heatmap.pdf", width = 24, height = 12)

pheatmap(af.mat, 
         cluster_cols = F,
         cluster_row = F,
         show_colnames = F,
         color=colorRampPalette(c("white","orange", "orangered", "red","darkred"))(100),
         annotation_col = chr,
         annotation_row = pop_basin, 
         annotation_colors = anno_colors,
         fontsize=15)

dev.off()


#difference in MAF between focal population and others
c=20
focal= "FOXR"
other <- order[! order %in% focal]

focal1="CXVL"
focal2="LTTV"
other <- order[! order %in% c(focal1, focal2)]

af %>% 
  filter(chr == c) %>%
  select(chr, snp, order) %>% 
  summary()


af %>% 
  filter(chr == c) %>% 
  select(chr, snp, order) %>% 
  mutate(
    focal= rowMeans(select(., !!!rlang::syms(focal))),
    other= rowMeans(select(., !!!rlang::syms(other)))) %>% 
  select(chr, snp, focal, other) %>% 
  summary()

af %>% 
  filter(chr == c) %>% 
  select(chr, snp, order) %>% 
  mutate(
    focal1=CXVL,
    focal2=LTTV,
    #focal= rowMeans(select(., !!!rlang::syms(focal))),
    other= rowMeans(select(., !!!rlang::syms(other)))) %>% 
  select(chr, snp, focal1, focal2, other) %>% 
  #select(chr, snp, focal, other) %>% 
  summary()










