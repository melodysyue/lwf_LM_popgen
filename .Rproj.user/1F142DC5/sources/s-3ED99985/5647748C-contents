library(tidyverse)
library(here)
library(ggsci)
library(gtools)


rm(list=ls())
##popmap
popmap <- read.table("./data/meta/lwf.N829.popmap.edited", sep='\t')
colnames(popmap) <- c("sample", "pop")
popmap$pop <- as.factor(popmap$pop)
order <- c( "FOXR", "RILY","LSRF", "PSTG","MENO", "WTFB","NMNL", "BBDN","MNST", 
            "NAUB", "EPFT", "CXVL", "LTTV", "EKRP", "GOOD", "IGPT", "MSKG")


popmap$pop=factor(popmap$pop, levels=order)

east <- c("NAUB", "EPFT", "CXVL", "LTTV", "EKRP", "GOOD", "IGPT", "MSKG")

popmap <- popmap %>% 
  filter(pop %in% east)


##ancestry proportion
names <- list.files(path="./data/admixture/", pattern = "lwf.east.pop.*.Q", full.names = F) 
files <- list.files(path="./data/admixture/", pattern = "lwf.east.pop.*.Q", full.names = T) %>% 
  lapply(read.table, header=F)

names <- str_remove_all(names, "lwf.east.pop.")
names <- str_remove_all(names, ".Q")
names <- paste0("K=",names)


fam <- read.table("./data/admixture/lwf.east.pop.fam",header=F) %>% 
  select(sample=V2)


admix=NULL #create an empty data frame to store data

for (i in 1: length(names)){
  k=files[[i]]
  npop <- names[[i]]
  k <- cbind(fam, k)
  k <- left_join(popmap, k, by="sample") %>% 
    mutate(k=npop) %>% 
    select(sample, pop, k, everything()) %>% 
    gather(key = "infered_clump", value = "prop_clump", 4:ncol(.))
  
  admix <- rbind(admix, k)
}

admix$k <- as.factor(admix$k)
admix$k <- factor(admix$k, levels=paste0("K=",seq(1,10,1)))

##plot admixture plot

pdf("./figures/figS2_admixture_east.pdf", width = 15, height = 9)
admix %>% 
  ggplot(aes(x = sample, y = as.numeric(prop_clump), fill = infered_clump))+
  geom_bar(stat = "identity", width=1.2) +
  facet_grid(k ~ pop, scales = "free_x", space = "free", switch = "both",
             labeller = labeller(k=k)) +
  theme_bw(base_size = 20) +
  scale_fill_simpsons()+
  theme(legend.position = "none", 
        axis.text = element_blank(), 
        axis.title = element_blank(),
        panel.spacing = unit(0,"lines"),
        panel.border = element_rect(fill=NA, size=0.5),
        strip.background = element_blank(),
        axis.line.x=element_blank(),
        axis.ticks=element_blank(),
        strip.text.x = element_text(angle=90),
        strip.text.y = element_text(angle=180))+
  scale_y_continuous(expand = c(0, 0))

dev.off()



