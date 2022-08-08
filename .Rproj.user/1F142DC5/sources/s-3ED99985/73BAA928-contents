library(tidyverse)
library(here)
library(ggsci)
library(gtools)
library(patchwork)


rm(list=ls())
##popmap
popmap <- read.table("./data/meta/lwf.N829.popmap.edited", sep='\t')
colnames(popmap) <- c("sample", "pop")
popmap$pop <- as.factor(popmap$pop)

order <- c( "FOXR", "RILY","LSRF", "PSTG","MENO", "WTFB","NMNL", "BBDN","MNST", 
            "NAUB", "EPFT", "CXVL", "LTTV", "EKRP", "GOOD", "IGPT", "MSKG")

popmap$pop=factor(popmap$pop, levels=order)

west <- c("FOXR", "RILY","LSRF", "PSTG","MENO", "WTFB","NMNL", "BBDN","MNST")
east <- c("NAUB", "EPFT", "CXVL", "LTTV", "EKRP", "GOOD", "IGPT", "MSKG")


#allpops
Q.all <- read.table("./data/admixture/lwf.N829.2.Q", header=F)
fam.all <- read.table("./data/admixture/lwf.N829.fam", header = F) %>% 
  select(sample=V2)
k.all <-cbind(fam.all, Q.all)
k.all <- left_join(popmap, k.all, by="sample") %>% 
  gather(key = "infered_clump", value = "prop_clump", 3:ncol(.))


p.all.ad <- k.all %>% 
  ggplot(aes(x = sample, y = as.numeric(prop_clump), fill = infered_clump))+
  geom_bar(stat = "identity", width=1.2) +
  facet_grid( ~ pop, scales = "free_x", space = "free", switch = "both") +
  theme_bw(base_size = 20) +
  scale_fill_simpsons()+
  labs(title="All Populations (K=2)")+
  theme(legend.position = "none",
        legend.title = element_blank(),
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

all.cv <- read.table("./data/admixture/choose.lwf.N829.txt",header=FALSE)
colnames(all.cv) <- c("K", "all.cv")
p.all.cv <- all.cv %>% 
  ggplot(aes(x=K, y=all.cv))+
  geom_line()+
  theme_classic(base_size = 20)+
  scale_x_continuous(breaks=seq(1,10,1))+
  ylab("CV error")




#west
Q.w <- read.table("./data/admixture/lwf.west.pop.2.Q", header=F)
fam.w <- read.table("./data/admixture/lwf.west.pop.fam", header = F) %>% 
  select(sample=V2)
k.w <-cbind(fam.w, Q.w)
popmap.w <- popmap %>% 
  filter(pop %in% west)
k.w <- left_join(popmap.w, k.w, by="sample") %>% 
  gather(key = "infered_clump", value = "prop_clump", 3:ncol(.))

p.w.ad <- k.w %>% 
  ggplot(aes(x = sample, y = as.numeric(prop_clump), fill = infered_clump))+
  geom_bar(stat = "identity", width=1.2) +
  facet_grid( ~ pop, scales = "free_x", space = "free", switch = "both") +
  theme_bw(base_size = 20) +
  scale_fill_simpsons()+
  labs(title="Northwestern Lake Michigan (K=2)")+
  theme(legend.position = "none",
        legend.title = element_blank(),
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

w.cv <- read.table("./data/admixture/chooseK.lwf.west.txt",header=FALSE)
colnames(w.cv) <- c("K", "w.cv")
p.w.cv <- w.cv %>% 
  ggplot(aes(x=K, y=w.cv))+
  geom_line()+
  theme_classic(base_size = 20)+
  scale_x_continuous(breaks=seq(1,10,1))+
  ylab("CV error")


#east
Q.e <- read.table("./data/admixture/lwf.east.pop.5.Q", header=F)
fam.e <- read.table("./data/admixture/lwf.east.pop.fam", header = F) %>% 
  select(sample=V2)
k.e <-cbind(fam.e, Q.e)
popmap.e <- popmap %>% 
  filter(pop %in% east)
k.e <- left_join(popmap.e, k.e, by="sample") %>% 
  gather(key = "infered_clump", value = "prop_clump", 3:ncol(.))

p.e.ad <- k.e %>% 
  ggplot(aes(x = sample, y = as.numeric(prop_clump), fill = infered_clump))+
  geom_bar(stat = "identity", width=1.2) +
  facet_grid( ~ pop, scales = "free_x", space = "free", switch = "both") +
  theme_bw(base_size = 20) +
  scale_fill_simpsons()+
  labs(title="Eastern Lake Michigan (K=5)")+
  theme(legend.position = "none",
        legend.title = element_blank(),
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

e.cv <- read.table("./data/admixture/chooseK.lwf.east.txt",header=FALSE)
colnames(e.cv) <- c("K", "e.cv")
p.e.cv <- e.cv %>% 
  ggplot(aes(x=K, y=e.cv))+
  geom_line()+
  theme_classic(base_size = 20)+
  scale_x_continuous(breaks=seq(1,10,1))+
  ylab("CV error")



pdf("./figures/figS1_admixture.pdf", width = 12, height = 9)

p.all.ad /p.w.ad /p.e.ad 

dev.off()


pdf("./figures/figS1_cv.pdf", width = 4, height = 9)

p.all.cv /p.w.cv /p.e.cv +
  plot_annotation(tag_levels = 'A')
dev.off()
