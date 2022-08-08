library(adegenet) #dapc
library(here)
library(tidyverse)
library(hierfstat)#basic.stats
library(poppr) #popsub, aboot
library(RColorBrewer)
library(ggtree)
library(ape) #write.tree
library(pegas)
library(zvau)
library(Matrix)
library(ggpubr)


rm(list=ls())

########################
## Make Genepop Object##
########################

#double check if your .gen file has a title, if it is blank, add a title. 
#Otherwise it will give you an error. Read.genepop file can take 10 min. 
#read.genepop reads .gen file and convert it into a genind object.

gen <- read.genepop("./data/islands/chr20_region_ol.gen", ncode=3L)
gen

s.keep <- indNames(gen) #to get the list of individuals

#add pop information
popmap <- read.table("./data/meta/lwf.N829.popmap.edited",header=F)
colnames(popmap)=c("sample","pop")
p.keep <- popmap[popmap$sample %in% s.keep,]
table(s.keep == p.keep$sample)
pop(gen)=p.keep$pop
gen

#double check
nPop(gen)
popNames(gen)
table(p.keep$pop)

#################
###   pwFst   ###
#################

##pairwise fst
ppfst <- genet.dist(gen, method="WC84")
ppfst <- as.data.frame(as.matrix(ppfst))
ppfst

ppfst.long <- ppfst %>% 
  as.data.frame() %>% 
  rownames_to_column("pop1") %>% 
  gather(pop2, pwfst, NMNL:RILY) %>% 
  drop_na() %>% 
  filter(pop1!=pop2)

##pairwise fst significance
hf <- genind2hierfstat(gen) 
pop.levels <- levels(gen$pop)
pop.levels

##for boot.ppfst to run properly, the pop column needs to be in numeric form and sorted.
pop <- sort(as.numeric(hf[,1])) 
pop
ppfst.ci<- boot.ppfst(data.frame(pop, hf[,-1]), nboot=1000, quant=c(0.025, 0.975))

pop.levels
ppfst.ci$ll
ppfst.ci$ul

colnames(ppfst.ci$ll) <- pop.levels
rownames(ppfst.ci$ll) <- pop.levels
colnames(ppfst.ci$ul) <- pop.levels
rownames(ppfst.ci$ul) <- pop.levels

#make symmetrical matrix 
ppfst.cill <- forceSymmetric(ppfst.ci$ll, uplo="U")
ppfst.ciul <- forceSymmetric(ppfst.ci$ul, uplo="U")


ppfst.cill <- ppfst.cill %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("pop1") %>% 
  gather(pop2, ci_ll, NMNL:RILY) %>% 
  drop_na()

ppfst.ciul <- ppfst.ciul %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("pop1") %>% 
  gather(pop2, ci_ul, NMNL:RILY) %>% 
  drop_na()

ppfst.p <- left_join(ppfst.cill, ppfst.ciul)


ppfst.summary <- left_join(ppfst.long, ppfst.p)

ppfst.summary <- ppfst.summary %>% 
  mutate(p = if_else(ci_ll<=0, "NS", "S"))


#pwfst heatmap
#to make the color scale consistent, we will use scale of (0, 0.4)
ppfst.summary$pwfst[ppfst.summary$pwfst<0]=0
summary(ppfst.summary)

#make sure lower limit is smaller than upper limit
table(ppfst.summary$ci_ll<=ppfst.summary$ci_ul)


write.table(ppfst.summary, "./data/islands/chr20_ol_pwfst.txt", quote=F, row.names = F)


