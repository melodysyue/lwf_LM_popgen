library(tidyverse)
library(ggpubr)
library(ggsci)
library(patchwork)

rm(list=ls())
##popmap
popmap <- read.table("./data/meta/lwf.N829.popmap.edited", sep='\t')
colnames(popmap) <- c("sample", "pop")
order <- c( "FOXR", "RILY","LSRF", "PSTG","MENO", "WTFB","NMNL", "BBDN","MNST", 
            "NAUB", "EPFT", "CXVL", "LTTV", "EKRP", "GOOD", "IGPT", "MSKG")

popmap$pop <- factor(popmap$pop, levels=order)

n=length(unique(popmap$pop))

popcol=setNames(pal_igv("default")(17), levels(popmap$pop))

west <- c("FOXR", "RILY","LSRF", "PSTG","MENO", "WTFB","NMNL", "BBDN","MNST")
east <- c("NAUB", "EPFT", "CXVL", "LTTV", "EKRP", "GOOD", "IGPT", "MSKG")



#all pops
  
t_all <- read.table("./data/plink_pca/lwf.eigenvec")
vals_all <- read.table("./data/plink_pca/lwf.eigenval")

colnames(t_all) <- c("pop", "sample", paste0("PC", seq(1,10,1)))
t_all$pop <- factor(t_all$pop, levels=order)
pops_all <- unique(t_all$pop)
col.subset_all <- popcol[names(popcol) %in% pops_all]

p_all <- t_all %>% 
  ggplot(aes (x =PC1, y = PC2, color=pop))+
  geom_point(size = 2, alpha =0.7)+
  scale_color_manual(values=col.subset_all)+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept=0, linetype="dashed")+
  theme_bw(base_size = 20)+
  theme(legend.position = "right",
        legend.title=element_blank())+
  labs(x = paste0("PC1 (", round(vals_all$V1[1],1), "%)"),
       y = paste0("PC2 (", round(vals_all$V1[2],1), "%)"),
       title=paste0("All Populations"))+
  guides(color=guide_legend(ncol=2, byrow=FALSE))



leg <- get_legend(p_all)
p_leg <- as_ggplot(leg)

#all pops by basin
t_all <- t_all %>% 
  mutate(basin=case_when(
    pop %in% west ~ "Northwest",
    pop %in% east ~ "East"))
t_all$basin <- factor(t_all$basin, levels=c("Northwest", "East"))

basincol <- setNames(c("steelblue", "orange"), levels(t_all$basin))




p_all_basin <- t_all %>% 
  ggplot(aes (x =PC1, y = PC2, color=basin))+
  geom_point(size = 2, alpha =0.7)+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept=0, linetype="dashed")+
  scale_color_manual(values=basincol)+ 
  theme_bw(base_size = 20)+
  theme(legend.position = c(0.3, 0.13),
        legend.title=element_blank(), 
        legend.background = element_rect(fill="white", color="black", size=0.5))+
  labs(x = paste0("PC1 (", round(vals_all$V1[1],1), "%)"),
       y = paste0("PC2 (", round(vals_all$V1[2],1), "%)"))+
  guides(color=guide_legend(nrow=1))



#west
t_w <- read.table("./data/plink_pca/lwf.west.eigenvec")
vals_w <- read.table("./data/plink_pca/lwf.west.eigenval")

colnames(t_w) <- c("pop", "sample", paste0("PC", seq(1,10,1)))
t_w$pop <- factor(t_w$pop, levels=order)
pops_w <- unique(t_w$pop)
col.subset_w <- popcol[names(popcol) %in% pops_w]

p_w <- t_w %>% 
  ggplot(aes (x =PC1, y = PC2, color=pop))+
  geom_point(size = 2, alpha =0.7)+
  scale_color_manual(values=col.subset_w)+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept=0, linetype="dashed")+
  theme_bw(base_size = 20)+
  labs(x = paste0("PC1 (", round(vals_w$V1[1],1), "%)"),
       y = paste0("PC2 (", round(vals_w$V1[2],1), "%)"),
       title="Northwest")+
  theme(legend.position = "none",
        legend.title=element_blank(),
        panel.border = element_rect(colour="steelblue", size=2),
        text=element_text(colour="steelblue"))

#east

t_e <- read.table("./data/plink_pca/lwf.east.eigenvec")
vals_e <- read.table("./data/plink_pca/lwf.east.eigenval")


colnames(t_e) <- c("pop", "sample", paste0("PC", seq(1,10,1)))
t_e$pop <- factor(t_e$pop, levels=order)
pops_e <- unique(t_e$pop)
col.subset_e <- popcol[names(popcol) %in% pops_e]

p_e <- t_e %>% 
  ggplot(aes (x =PC1, y = PC2, color=pop))+
  geom_point(size = 2, alpha =0.7)+
  scale_color_manual(values=col.subset_e)+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept=0, linetype="dashed")+
  theme_bw(base_size = 20)+
  labs(x = paste0("PC1 (", round(vals_e$V1[1],1), "%)"),
       y = paste0("PC2 (", round(vals_e$V1[2],1), "%)"),
       title="East")+
  theme(legend.position = "none",
        legend.title=element_blank(),
        panel.border = element_rect(colour="orange", size=2),
        text=element_text(colour="orange"))




#northeastern

t_ne <- read.table("./data/plink_pca/lwf.ne.eigenvec")
vals_ne <- read.table("./data/plink_pca/lwf.ne.eigenval")

colnames(t_ne) <- c("pop", "sample", paste0("PC", seq(1,10,1)))
t_ne$pop <- factor(t_ne$pop, levels=order)
pops_ne <- unique(t_ne$pop)
col.subset_ne <- popcol[names(popcol) %in% pops_ne]

p_ne <- t_ne %>% 
  ggplot(aes (x =PC1, y = PC2, color=pop))+
  geom_point(size = 2, alpha =0.7)+
  scale_color_manual(values=col.subset_ne)+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept=0, linetype="dashed")+
  theme_bw(base_size = 20)+
  labs(x = paste0("PC1 (", round(vals_ne$V1[1],1), "%)"),
       y = paste0("PC2 (", round(vals_ne$V1[2],1), "%)"),
       title="Northeast")+
  theme(legend.position = "none",
        legend.title=element_blank(),
        panel.border = element_rect(colour="salmon", size=2),
        text=element_text(colour="salmon"))




#save the plot
both <- ((p_all + theme(legend.position = "none")) + p_all_basin) /
  (p_w + p_e + p_ne) 


pdf("./figures/fig1BC_pca.pdf", width = 12, height = 9)
both
dev.off()


pdf("./figures/fig1_legend.pdf", width = 6, height = 9)
p_leg
dev.off()





