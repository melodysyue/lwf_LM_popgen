library(tidyverse)

rm(list=ls())

all <- read.table("./data/neutral_outlier/lwf.pcadapt.allpops.qvalues.txt", header=T)
west <-read.table("./data/neutral_outlier/lwf.pcadapt.west.qvalues.txt", header=T)
east <- read.table("./data/neutral_outlier/lwf.pcadapt.east.qvalues.txt", header=T)
ne <- read.table("./data/neutral_outlier/lwf.pcadapt.ne.qvalues.txt", header=T)

all$data <- "All populations"
west$data <- "Northwest"
east$data <- "East"
ne$data <- "Northeast"

dd <- rbind(all, west, east, ne)
dd$chr <- as.factor(dd$chr)
dd$chr <- factor(dd$chr, levels=seq(1, length(unique(dd$chr))))

dd$data <- as.factor(dd$data)
dd$data <- factor(dd$data, levels=c("All populations", "Northwest", "East", "Northeast"))


dd <- dd %>% 
  mutate(ol=qval<0.01) %>% 
  filter(!chr==40)

axis.set <- dd %>% 
  group_by(chr) %>% 
  summarize(center = (max(chr.mbp) + min(chr.mbp)) / 2) %>% 
  filter(!chr==40)



#islands
dd.r <- read.table("./data/islands/sigRegions_start_end.txt", header=T)
dd.r
duptimes <- c(1,1,2,1,1,1)
idx <- rep(1:nrow(dd.r), duptimes)  
dd.r <- dd.r[idx,]
dd.r

dd.r$data <- c("All populations", "East", "East", "Northeast", "All populations", "East", "Northwest")
dd.r$data <- factor(dd.r$data, levels=c("All populations", "Northwest", "East", "Northeast"))
dd.r$chr <- as.factor(dd.r$chr)

  
dd.r <- left_join(dd.r, dd, by=c("chr", "start"="pos", "data")) %>% 
  select(data, chr, start=chr.mbp, end)
dd.r <- left_join(dd.r, dd, by=c("chr", "end"="pos", "data")) %>% 
  select(data, chr, start, end=chr.mbp) %>% 
  unique() %>% 
  mutate(size=end-start)



p.qval <- ggplot(dd) +
  geom_point(aes(x = chr.mbp, y = logq_plot, color = chr), size=2) +
  scale_x_continuous(label = axis.set$chr, breaks = axis.set$center) +
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0, 10))+
  scale_color_manual(values = rep(c("dimgray","darkgrey"), 20)) +
  geom_hline(yintercept = -log10(0.01), linetype="dashed", color="red", size=2)+
  geom_rect(data=dd.r, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill="purple", alpha=0.6)+
  geom_point(data=subset(dd, ol=="TRUE"), aes(x=chr.mbp, y=logq_plot), col="orange", size=2)+
  ylab("-log10(q-value)")+
  theme_bw(base_size = 25)+
  facet_wrap(~data, dir="v", scales="free", ncol=1, strip.position="right")+
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(vjust = 0.5, size=15),
    strip.background = element_blank(),
    axis.title.x = element_blank()
  )



ggsave(filename = "./figures/fig3.manhattan.png", p.qval,  width=24, height=16, units="in")

pdf("./figures/fig3.manhattan.pdf", width = 24, height = 16)
p.qval
dev.off()

