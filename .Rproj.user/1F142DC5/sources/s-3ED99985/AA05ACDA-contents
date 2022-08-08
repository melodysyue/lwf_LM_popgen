library(tidyverse)

rm(list=ls())


#pwfst 
mt <- read.table("./data/fst/pwfst.txt", header=T)

#value range by basin
west <- c("FOXR", "RILY","LSRF", "PSTG","MENO", "WTFB","NMNL", "BBDN","MNST")
east <- c("NAUB", "EPFT", "CXVL", "LTTV", "EKRP", "GOOD", "IGPT", "MSKG")

mt %>% 
  filter( var1 %in% west & var2 %in% west) %>% 
  filter(var1 != var2) %>% 
  summary()

mt %>% 
  filter( var1 %in% east & var2 %in% east) %>% 
  filter(var1 != var2) %>% 
  view()



bc_a <- 0.05/136 #multi-testing correction
max <- round(max(mt$fst),3)

#order
order <- c( "FOXR", "RILY","LSRF", "PSTG","MENO", "WTFB","NMNL", "BBDN","MNST", 
            "NAUB", "EPFT", "CXVL", "LTTV", "EKRP", "GOOD", "IGPT", "MSKG")
mt$var1 <- factor(mt$var1, levels=order)
mt$var2 <- factor(mt$var2, levels=order)

text.col <- c(rep("steelblue", length(west)), rep("orange", length(east)))


pdf("./figures/figS3_pwfst.pdf", width = 12, height = 9)
mt%>%
  ggplot() +
  geom_tile(aes(x = var1, y = var2, fill = fst), color = "white")+
  geom_point(data=subset(mt, pval>bc_a), aes(x = var1, y = var2), shape=4, size=6)+
  scale_fill_gradient(low = "white", high = "red", 
                      limit = c(0, max), 
                      breaks=c(0, 0.005, 0.010, max),
                      name=expression(italic("F")[ST])) +
  theme_classic(base_size = 20)+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
        axis.text = element_text(colour=text.col),
        axis.title = element_blank(),
        axis.ticks = element_blank())+
  coord_fixed()

dev.off()


