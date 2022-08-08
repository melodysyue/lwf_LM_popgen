library(tidyverse)
library(vcfR)
library(pheatmap)
library(ggsci)
library(scales)
library(patchwork)
library(ggplotify)

rm(list=ls())

popmap <- read.table("./data/meta/lwf.N829.popmap.edited", header=F)
colnames(popmap) <- c("sample", "pop")
head(popmap)

order <- c( "FOXR", "RILY","LSRF", "PSTG","MENO", "WTFB","NMNL", "BBDN","MNST", 
            "NAUB", "EPFT", "CXVL", "LTTV", "EKRP", "GOOD", "IGPT", "MSKG")


popmap$pop <- factor(popmap$pop, levels=order)
popcol=setNames(pal_igv("default")(17), levels(popmap$pop))

vcf <- read.vcfR("./data/vcf/chr20_region_ol.recode.vcf")
c <- "20"


##################
###lostruct MDS###
##################
region <- read.table("./data/islands/sigRegions_start_end.txt", header=T)
mds.coords.loci.win <- read.table("./data/lostruct/20/20_lostruct_output.txt", header=T)

start <- region %>% 
  filter(chr==c) %>% 
  pull(start)

end <- region %>% 
  filter(chr==c) %>% 
  pull(end)

region.win <- mds.coords.loci.win %>% 
  select(window, rowID, chr, pos) %>% 
  filter(pos>=start & pos<=end) %>% 
  pull(window) %>% 
  unique()

##MDS plot for MDS1
i <- 1 #MDS1
n.sd <- 3
cutoff <- mean(mds.coords.loci.win[,i])+c(-1,1)*n.sd*sd(mds.coords.loci.win[,i]) 
p_lostruct <- mds.coords.loci.win %>% 
  ggplot(aes(x=window, y= mds.coords.loci.win[,i])) + 
  geom_point()+
  geom_hline(yintercept = cutoff, linetype="dashed", color="red", size=0.5) +
  annotate("rect", fill="purple", alpha=0.3, xmin=min(region.win), xmax=max(region.win), ymin=-Inf, ymax=Inf)+
  theme_bw(base_size = 20)+
  xlab(paste0("Window index on Chromosome ", chr))+
  ylab(paste0("MDS",i))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())  

p_lostruct


##############
####genind####
##############
gen <- vcfR2genind(vcf)
s.keep <- indNames(gen)
popmap.sorted <- popmap[match(s.keep, popmap$sample),]
table(s.keep==popmap.sorted$sample)
pop(gen) <- popmap.sorted$pop
gen

############
####DAPC####
############

# cluster/dapc
set.seed(29)
grp <- find.clusters(gen, n.pc=1, n.clust = 2)
cluster <- grp$grp %>% 
  as.data.frame() %>% 
  rownames_to_column("sample")
colnames(cluster)=c("sample", "cluster")

table(cluster$cluster)
cluster <- cluster %>% 
  mutate(grp=case_when(
    cluster==1 ~ "0",
    cluster==2 ~ "1",
  ))

cluster$grp <- as.factor(cluster$grp)
cluster$grp <- factor(cluster$grp, levels=c("0", "1"))
table(cluster$grp)
table(cluster$cluster)

cluster.pop <- left_join(cluster, popmap, by="sample")
head(cluster.pop)
cluster.pop %>% 
  filter(grp==1) %>% 
  group_by(pop) %>% 
  summarize(n=n())

table(popmap$pop)
###################
####PCA+cluster####
###################

trans <- scaleGen(gen, NA.method="mean")
pca_obj <- dudi.pca(trans,cent=FALSE,scale=FALSE,scannf=FALSE, nf=3)
eigenvects_to_plot=pca_obj$li
variance_explained=(pca_obj$eig)/sum(pca_obj$eig)
xlabel=paste0("PC",1," (",round(variance_explained[1],3)*100,"%)")
ylabel=paste0("PC",2," (",round(variance_explained[2],3)*100,"%)")
eigenvects_to_plot <- eigenvects_to_plot %>% 
  rownames_to_column("sample")

eigenvects_to_plot <- left_join(eigenvects_to_plot, popmap, by="sample")
eigenvects_to_plot <- left_join(eigenvects_to_plot, cluster, by="sample")
head(eigenvects_to_plot)
eigenvects_to_plot$pop <- factor(eigenvects_to_plot$pop, levels=order)



pca_pop <- 
  ggplot(eigenvects_to_plot, aes(x=Axis1, y=Axis2, color=pop)) + #don't use data$column inside aes
  geom_point(size=3, alpha=0.8)+
  labs(x=xlabel, y=ylabel,
       col="Populations")+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept=0, linetype="dashed")+
  theme_bw(base_size=20)+
  scale_color_manual(values=popcol)+
  theme(legend.position = "bottom")+
  guides(colour=guide_legend(nrow=2)) #legend as one row

pca_pop

pca_cluster <- 
  ggplot(eigenvects_to_plot, aes(x=Axis1, y=Axis2, color=grp)) + #don't use data$column inside aes
  geom_point(size=3, alpha=0.8)+
  labs(x=xlabel, y=ylabel,
       col="Cluster")+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept=0, linetype="dashed")+
  theme_bw(base_size=20)+
  scale_color_manual(values=c("red", "blue"))+
  theme(legend.position = "bottom")+
  guides(colour=guide_legend(nrow=1)) #legend as one row

pca_cluster

########################
####genotype heatmap####
########################

#extract genotypes in numeric format
genotypes_numeric<-extract.gt(vcf, element = "GT", mask = FALSE, as.numeric = FALSE,
                              return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE,
                              convertNA = TRUE)

head(genotypes_numeric) #take a peak and the genotype is "|"

dim(genotypes_numeric)
table(is.na(genotypes_numeric))


#function to convert genotypes to 0, 1, 2 
#0 (homozygote for ref allele, 0|0)
#1 (hetero, 0|1)
#2 (homozygote for alt allele, 1|1) 
convertGenos<-function(genotype){
  alleles<-as.numeric(str_split_fixed(genotype,"\\/",2))
  genoCode<-alleles[1]+alleles[2]
  return(genoCode)
}

#convert to 0,1,2 format
genosCode<-apply(genotypes_numeric,1:2,convertGenos)
genosCode<-t(genosCode)

#input format for heatmap
heatmap<-genosCode
heatmap[is.na(heatmap)] <- -1 #replace missing info with -1
heatmap<-apply(heatmap,1:2,function(x) as.numeric(x))
dim(heatmap)
head(heatmap)[,1:6]

#manipulate heatmap matrix
heatmap <- heatmap %>% 
  as.data.frame() %>% 
  rownames_to_column("sample")
head(heatmap)[,1:6]

rownames(heatmap) <- heatmap$sample
heatmap <- heatmap[,-1]


#row annotation (for samples) and order row by cluster
popmap.cluster <- left_join(popmap.sorted, cluster, by="sample")

row.order <- popmap.cluster %>% 
  arrange(grp,pop) %>% 
  pull(sample)

heatmap <- heatmap %>% 
  as.data.frame() %>% 
  rownames_to_column("sample")
heatmap.sorted <- heatmap[match(row.order, heatmap$sample),]
table(heatmap.sorted$sample==row.order)
rownames(heatmap.sorted) <- heatmap.sorted$sample
heatmap.sorted <- heatmap.sorted[,-1]

row_anno <- popmap.cluster %>% 
  select(sample, pop, cluster=grp)
row_anno_sorted <- row_anno[match(row.order, row_anno$sample),]
table(row_anno_sorted$sample==row.order)
rownames(row_anno_sorted) <- row_anno_sorted$sample
row_anno_sorted <- row_anno_sorted[,-1]
row_anno_sorted$cluster <- as.factor(row_anno_sorted$cluster)
row_anno_sorted$pop <- factor(row_anno_sorted$pop, levels=order)

levels(row_anno_sorted$cluster)



###color
cluster_colors <- list(cluster=c("red", "blue"))
names(cluster_colors$cluster) <- c("0","1")


pop_colors <- list(pop=popcol)
names(pop_colors$pop) <- order

mycol <- append(cluster_colors,pop_colors)


#plot it
p_heatmap <- pheatmap(heatmap.sorted, 
                      annotation_row = row_anno_sorted,
                      annotation_colors = mycol,
                      cluster_rows = FALSE,
                      cluster_cols = FALSE,
                      show_colnames = FALSE,
                      show_rownames = FALSE,
                      color=c("lightblue","steelblue", "darkblue"),
                      legend_breaks = c(0,1,2),
                      legend_labels =c("REF homozygous","heterozygous","ALT homozygous"),
                      fontsize = 15)


p_heatmap <- as.ggplot(p_heatmap)


###put them all together
p_pca <- ((pca_pop + theme(legend.position = "none")) | (pca_cluster + theme(legend.position = "none")))


pdf(paste0("./figures/fig5_chr", chr, "_lostruct_pca_genoheatmap.pdf"), width = 20, height = 18)
p_lostruct/ p_pca / p_heatmap + 
  plot_layout(ncol=1, heights = c(1,2,4))+
  plot_annotation(tag_levels = "A")
dev.off()