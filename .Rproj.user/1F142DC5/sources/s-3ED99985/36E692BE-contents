library(tidyverse)
library(pcadapt)
library(qvalue)


rm(list=ls())

#input
popmap <- read.table("./data/plink/lwf.ne.pop.fam", header=F)
popmap <- popmap %>% 
  select(sample=V2, pop=V1)
head(popmap)


pd <- read.pcadapt("./data/plink/lwf.ne.pop.bed", type="bed")
pd.thin <- read.pcadapt("./data/plink/lwf.pcadapt.ne.thinned.bed", type="bed")

snps <- read.table("./data/plink/lwf.ne.pop.bim", header = F)


#determine K
x <- pcadapt(input=pd.thin, K=20)

pdf("./figures/lwf.pcadapt.ne.thin.K.pdf")
plot(x, option = "screeplot")
plot(x, option = "scores", pop=popmap$pop) 
plot(x, option = "scores", pop=popmap$pop, i=3, j=4) 
dev.off()


#pcadapt
y <- pcadapt(input=pd, K=2, min.maf = 0.05)
qval.cutoff <- 0.01

qval <- qvalue(y$pvalues)$qvalues
pcadapt.ol <- which(qval<qval.cutoff)

pcadapt.ol.snps <- snps %>% 
  rownames_to_column("rowID") %>% 
  filter(rowID %in% pcadapt.ol) %>% 
  pull(V2)

length(pcadapt.ol.snps)

write.table(pcadapt.ol.snps, paste0("./data/neutral_outlier/lwf.pcadapt.ol.qval", qval.cutoff, ".ne.snps"),quote=F, row.names = F)

#spread across how many chromosomes
snps %>% 
  filter(V2 %in% pcadapt.ol.snps) %>% 
  group_by(V1) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n)) %>% 
  view()


#q values
snps.pval <- cbind(snps, qval=qval)
colnames(snps.pval) <- c("chr", "snp", "dummy", "pos", "minor", "major", "qval")
head(snps.pval)


#calculate cumulative pos
snps.pval$cumpos <- NA 
s <- 0
nbp <- c()

for (i in unique(snps.pval$chr)){
  nbp[i] <- max(snps.pval[snps.pval$chr == i,]$pos)
  snps.pval[snps.pval$chr == i,"cumpos"] <- snps.pval[snps.pval$chr == i,"pos"] + s
  s <- s + nbp[i]
}


snps.pval <- snps.pval %>% 
  mutate(chr.mbp=cumpos/1e6,
         logq=-log10(qval)) %>% 
  mutate(logq_plot=if_else(is.na(logq), 0.001, logq)) #replace NA with a very small value for plotting


write.table(snps.pval, "./data/neutral_outlier/lwf.pcadapt.ne.qvalues.txt", quote=F, row.names = F)



