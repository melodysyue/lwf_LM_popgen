library(lostruct)
library(tidyverse)

rm(list=ls())
#################
###Instruction###
#################
#lostruct.r takes three output file, *.geno, *.popmap, *.loci. 
#The script includes two sections
#Input and parameters section: change for your data;
#Lostruct section: simply run this chunk of codes, output will be saved in the corresponding chromosome folder;
#Output: *.pdf: MDS plot for each MDS axis;
#output: *_lostruct_output.txt: mds coordinates;
#output: *_lostruct_input_matrix.txt;


##########################
###Input and parameters###
##########################

chr <- 20
##parameters 
n.sd <- 3 #outlier window cutoff is 3 sd;
ws=50 #window size of 50 SNPs.The last window will be omitted as it is less than window size 
mds=40 #40 MDS axes.


##create a folder for each chromosome
if(!file.exists(paste0("./data/lostruct/",chr))){
  dir.create(paste0("./data/lostruct/",chr))
}

##input format: numeric genotype matrix. 
gen <- read.table(paste0("./data/lostruct/",chr,".geno"), header=FALSE)

##add individual ID as row names
popmap <- read.table("./data/meta/lwf.N829.popmap.edited", header=FALSE)
colnames(popmap) <- c("sample", "pop")
popmap$ind <- paste0(popmap$pop,"_", popmap$sample)
rownames(gen) <- popmap$ind

##add SNP pos as col names
loci <- read.table(paste0("./data/lostruct/",chr, ".loci"),header=FALSE)
colnames(loci) <- c("chr","pos")
colnames(gen) <- loci$pos

##replace NA with the most frequent genotype
gen <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))


##transpose it so rows as SNPs, columns as individuals.
snp <- t(gen)
mat <- as.matrix(snp)

##############
###lostruct###
##############

##run core algorithm of lostruct
pcs <- eigen_windows(mat, k=2, win=ws) 
pcdist <- pc_dist(pcs,npc=2)
mds <- cmdscale(pcdist, eig=TRUE, k=mds)
mds.coords <- mds$points
colnames(mds.coords) <- paste0("MDS", 1:ncol(mds.coords))
mds.coords <- mds.coords %>% 
  as.data.frame() %>% 
  mutate(window=seq(1,nrow(mds.coords),1))


##associate window index with SNP ID.
loci.win <- loci %>%
  rownames_to_column("rowID") %>% 
  mutate(rowID=as.numeric(rowID)) %>% 
  mutate(window=floor((rowID-1)/ws + 1))
mds.coords.loci.win <- left_join(mds.coords, loci.win, by="window")

#associate with candiate region
region <- read.table("./data/islands/sigRegions_start_end.txt", header=T)
start <- region %>% 
  filter(chromosome==chr) %>% 
  pull(start)

end <- region %>% 
  filter(chromosome==chr) %>% 
  pull(end)

region.win <- mds.coords.loci.win %>% 
  select(window, rowID, chr, pos) %>% 
  filter(pos>=start & pos<=end) %>% 
  pull(window) %>% 
  unique()

##MDS plot
pdf(paste0("./figures/lostruct_chr",chr,"_ws",ws,"_sd", n.sd, "_mds40.pdf"), width = 12, height = 4)
for (i in 1:40){
  #outlier window cutoff
  cutoff <- mean(mds.coords[,i])+c(-1,1)*n.sd*sd(mds.coords[,i]) 
  p <- mds.coords %>% 
    ggplot(aes(x=window, y= mds.coords[,i])) + 
    geom_point()+
    geom_hline(yintercept = cutoff, linetype="dashed", color="red", size=0.5) +
    annotate("rect", fill="purple", alpha=0.3, xmin=min(region.win), xmax=max(region.win), ymin=-Inf, ymax=Inf)+
    theme_bw(base_size = 20)+
    xlab("Window index on chromosomes")+
    ylab(paste0("MDS",i))+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())  
  print(p)
}
dev.off()

##save output
write.table(mds.coords.loci.win, paste0("./data/lostruct/", chr, "/",chr,"_lostruct_output.txt"),quote=F, row.names = F)
write.table(mat, paste0("./data/lostruct/",chr, "/",chr,"_lostruct_input_matrix.txt"), quote=F, row.names = T)
  


