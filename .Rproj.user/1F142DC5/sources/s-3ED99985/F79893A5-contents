library(vegan)
library(tidyverse)
library(reshape2)
library(car)
library(scales)
library(usedist)
library(patchwork)
library(ggpubr)


rm(list=ls())


##Load the data
wdist <- read.table("./data/meta/lwf_waterDist.txt", header=T)
pwfst<- read.table("./data/fst/pwfst.txt", header=T)


###mantel test###
# Prepare a data for the Mantel Test
# mantel function unfortunately only accepts a matrix as input
wdist.d <- as.dist(wdist)


pwfst.d <- pwfst %>% 
  mutate(linearized_fst = (fst/(1-fst))) %>% 
  dplyr::select(var1,var2, linearized_fst) %>% 
  spread(var2, linearized_fst) %>% 
  column_to_rownames("var1") %>% 
  as.matrix()


#make sure two matrixes have the same order
pwfst.d <- pwfst.d[labels(wdist.d), labels(wdist.d)]
pwfst.d <- as.dist(pwfst.d)

labels(wdist.d)
labels(pwfst.d)

table(labels(wdist.d)==labels(pwfst.d))

#mantel test
set.seed(29)
mantel(wdist.d, pwfst.d, method="pearson", permutations=10000)


#summary(lm(pwfst.d ~ wdist.d))

###plotting###
#put tow datasets together
d1 <- melt(as.matrix(wdist.d), varnames=c("row","col")) %>% 
  rename(water.dist=value)
d2 <- melt(as.matrix(pwfst.d), varnames=c("row","col")) %>% 
  rename(linearfst=value)


dd <- left_join(d1,d2,by=c("row", "col"))


west <- c("FOXR", "RILY","LSRF","PSTG", "WTFB", "MENO","NMNL","BBDN","MNST")
east <- c("NAUB","EPFT", "CXVL", "LTTV","IGPT","GOOD","EKRP","MSKG")

dd <- dd %>% 
  mutate(type=case_when(
    row %in% west & col %in% west ~ "Among northwestern populations",
    row %in% east & col %in% east ~ "Among eastern populations", 
    TRUE ~ "Between northwestern and eastern populations"
  ))

dd$type <- factor(dd$type, levels=c("Among northwestern populations", "Between northwestern and eastern populations", "Among eastern populations"))


dd <- dd %>% 
  mutate(pair=paste0(row,"_",col)) %>% 
  filter(row!=col) %>% #remove self pairs
  column_to_rownames("pair") %>% 
  dplyr::select(-row, -col) %>% 
  unique()  #remove duplicate pairs;


#mantel test by basin
table(labels(wdist.d)==labels(pwfst.d))

#west basin
w.wdist <- dist_subset(wdist.d, west)
w.pwfst <- dist_subset(pwfst.d, west)
table(labels(w.wdist)==labels(w.pwfst))
mantel(w.wdist, w.pwfst, method="pearson", permutations=10000)


#east basin
e.wdist <- dist_subset(wdist.d, east)
e.pwfst <- dist_subset(pwfst.d, east)
table(labels(e.wdist)==labels(e.pwfst))
mantel(e.wdist, e.pwfst, method="pearson", permutations=10000)


#mantle test without MSKG
noMSKG <- c("FOXR", "RILY","LSRF","PSTG", "WTFB", "MENO","NMNL","BBDN","MNST",
                    "NAUB","EPFT", "CXVL", "LTTV","IGPT","GOOD","EKRP")

noMSKG.wdist <- dist_subset(wdist.d, noMSKG)
noMSKG.pwfst <- dist_subset(pwfst.d, noMSKG)
table(labels(noMSKG.wdist)==labels(noMSKG.pwfst))
mantel(noMSKG.wdist, noMSKG.pwfst, method="pearson", permutations=10000)


#overall IBD
p1 <- dd %>% 
  ggplot(aes(x=water.dist, y=linearfst)) + #specify dataframe
  geom_point() +
  scale_x_continuous(labels=comma)+
  geom_smooth(method="lm", se=FALSE, color="red")+
  ylab(expression(italic(F[ST]/(1-F[ST])))) +   
  xlab("In-Water Distance \n (calculated as least-cost path distance)")  +
  ggtitle("All Populations")+
  annotate("text",
           label=c("paste(italic(r), \" = 0.468\")",
                   "paste(italic(p), \" = 0.005\")"),
           parse=TRUE,
           x=c(-Inf,-Inf), y=c(0.0123,0.0115), size=8, hjust=-0.5, color="red") +
  theme_bw(base_size = 20)


#IBD by basin
p2 <- dd %>% 
  ggplot(aes(x=water.dist, y=linearfst, col=type)) + #specify dataframe
  geom_point() +
  scale_color_manual(values=c("steelblue", "seagreen", "orange"))+
  scale_x_continuous(labels=comma)+
  geom_smooth(method="lm", se=FALSE)+
  labs(col="Comparison")+
  ggtitle("By Lake Region (MSKG included)")+
  theme_bw(base_size = 20)+
  theme(axis.title = element_blank(),
        legend.position = "bottom")+
  guides(color=guide_legend(ncol=1))



#IBD without MSKG
p3 <- dd %>% 
  rownames_to_column("pair") %>% 
  filter(!str_detect(pair, "MSKG")) %>% 
  ggplot(aes(x=water.dist, y=linearfst, col=type)) +
  geom_point() +
  scale_color_manual(values=c("steelblue", "seagreen", "orange"))+
  scale_x_continuous(labels=comma)+
  geom_smooth(method="lm", se=FALSE)+
  labs(col="Comparison")+
  ggtitle("By Lake Region (MSKG removed)")+
  theme_bw(base_size = 20)+
  theme(axis.title = element_blank(),
        legend.position = "bottom")+
  guides(color=guide_legend(ncol=1))


pdf("./figures/fig2_lwf_ibd.pdf", width = 16, height = 12)

(p1 | (p2 / p3 / guide_area() + plot_layout(heights = c(4,4,1), guides='collect'))) + 
  plot_annotation(tag_levels = "A")

dev.off()


