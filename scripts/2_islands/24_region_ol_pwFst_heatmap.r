library(tidyverse)

rm(list=ls())


#input
list <- list.files(path="./data/islands/", pattern="*_ol_pwfst.txt", full.names=TRUE, recursive = TRUE) %>% 
  sort()

chrs <- str_remove_all(list, "./data/islands/") %>% 
  sort()
chrs <- str_remove_all(chrs, "_ol_pwfst.txt")


out <- list %>%
  map_df(read.table, header=T) 

head(out)

out$chr <- rep(chrs, each=nrow(out)/6)

summary(out) #use scale (0, 0.5)


#order
order <- c( "FOXR", "RILY","LSRF", "PSTG","MENO", "WTFB","NMNL", "BBDN","MNST", 
            "NAUB", "EPFT", "CXVL", "LTTV", "EKRP", "GOOD", "IGPT", "MSKG")

out$pop1 <- factor(out$pop1, levels=order)
out$pop2 <- factor(out$pop2, levels=order)

out$chr <- as.factor(out$chr)
out$chr <- factor(out$chr, 
                  levels=c("chr4", "chr7", "chr10", "chr11", "chr18", "chr20"),
                  labels=c("Chromosome 4", "Chromosome 7", "Chromosome 10", "Chromosome 11", "Chromosome 18", "Chromosome 20"))


west <- c("FOXR", "RILY","LSRF", "PSTG","MENO", "WTFB","NMNL", "BBDN","MNST")
east <- c("NAUB", "EPFT", "CXVL", "LTTV", "EKRP", "GOOD", "IGPT", "MSKG")
text.col <- c(rep("steelblue", length(west)), rep("orange", length(east)))


pdf("./figures/fig4_region_ol_pwfst.pdf", width = 18, height = 12)

out%>%
  ggplot() +
  geom_tile(aes(x = pop1, y = pop2, fill = pwfst), color = "white")+
  geom_point(data=subset(out, p=="NS"), aes(x = pop1, y = pop2), shape=4, size=4)+
  scale_fill_gradient(low = "white", high = "red", 
                      limit = c(0, 0.5), 
                      breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
                      name=expression(italic("F")[ST])) +
  theme_classic(base_size = 25)+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, vjust=0.5),
        axis.text = element_text(colour=text.col),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust=0.5),
  )+
  coord_fixed()+
  facet_wrap(~chr)

dev.off()


