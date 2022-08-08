library(tidyverse)
library(ggpubr)

rm(list=ls())

###pop order
west <- c("FOXR", "RILY","LSRF", "PSTG","MENO", "WTFB","NMNL", "BBDN","MNST")
east <- c("NAUB", "EPFT", "CXVL", "LTTV", "EKRP", "GOOD", "IGPT", "MSKG")
order <- c( "FOXR", "RILY","LSRF", "PSTG","MENO", "WTFB","NMNL", "BBDN","MNST", 
            "NAUB", "EPFT", "CXVL", "LTTV", "EKRP", "GOOD", "IGPT", "MSKG")

text.col <- c(rep("steelblue", length(west)), rep("orange", length(east)))


###full pwfst
full <- read.table("./data/fst/pwfst.txt", header=T)
full$var1 <- factor(full$var1, levels=order)
full$var2 <- factor(full$var2, levels=order)
full <- full %>% 
  filter(var1 != var2)


###islands pwfst
list <- list.files(path="./data/islands/", pattern="*_ol_pwfst.txt", full.names=TRUE, recursive = TRUE) %>% 
  sort()

chrs <- str_remove_all(list, "./data/islands/") %>% 
  sort()
chrs <- str_remove_all(chrs, "_ol_pwfst.txt")

out <- list %>%
  map_df(read.table, header=T) 

out$chr <- rep(chrs, each=nrow(out)/6)


out$pop1 <- factor(out$pop1, levels=order)
out$pop2 <- factor(out$pop2, levels=order)
out$chr <- as.factor(out$chr)
out$chr <- factor(out$chr, 
                  levels=c("chr4", "chr7", "chr10", "chr11", "chr18", "chr20"),
                  labels=c("Chromosome 4", "Chromosome 7", "Chromosome 10", "Chromosome 11", "Chromosome 18", "Chromosome 20"))



###compare full vs islands

head(full)
head(out)

full$pairs <- paste0(full$var1, "_", full$var2)
out$pairs <- paste0(out$pop1, "_", out$pop2)

table(unique(full$pairs)==unique(out$pairs))



full_pwfst <- full %>% 
  select(pairs, full=fst)


island_pwfst <- out %>% 
  select(pairs, islands=pwfst, chr)

head(full_pwfst)
head(island_pwfst)

both <- left_join(island_pwfst, full_pwfst, by="pairs") %>% 
  gather(source, value, -pairs, -chr)


both$source <- as.factor(both$source)
both$source <- factor(both$source, levels=c("islands", "full"))


p <- both %>% 
  ggplot(aes(x=source, y=value, fill=source)) +
  geom_boxplot(alpha=0.5)+
  scale_fill_manual(values=c("red", "dimgray"), labels=c("Pcadapt outliers", "All SNPs"))+
  stat_compare_means(label.x=1.8,
                     label.y=0.45)+
  ylab(expression(Pairwise~italic(F[ST])))+
  theme_classic(base_size = 20)+
  theme(axis.title.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=1),
        legend.position = "bottom",
        legend.title = element_blank())+
  facet_wrap(~chr)

pdf("./figures/figS4_pwfst_islands_vs_full.pdf", width = 12, height = 9)
p
dev.off()





