library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(tools)
library(ggrepel)
library(ggspatial)
library(maps)
library(ggsci)
library(grid)
library(s2)

rm(list=ls())

sf_use_s2(FALSE) #switch off S2, otherwise it won't work


world <- ne_countries(scale="medium", returnclass = "sf")
class(world)
loc <- read.csv("./data/meta/collections_lat_longs_afterfiltering_corrected.csv",header=TRUE)

order <- c( "FOXR", "RILY","LSRF", "PSTG","MENO", "WTFB","NMNL", "BBDN","MNST", 
            "NAUB", "EPFT", "CXVL", "LTTV", "EKRP", "GOOD", "IGPT", "MSKG")

loc <- loc %>% 
  mutate(pop=code)
loc$pop <- factor(loc$pop, levels=order)


states <-st_as_sf(maps::map("state", plot=FALSE, fill=TRUE))
states <- cbind(states, st_coordinates(st_centroid(states)))

states$ID <- str_replace(states$ID, "michigan", "MICHIGAN")
states$ID <- str_replace(states$ID, "wisconsin", "WISCONSIN")
states$ID <- str_replace(states$ID, "illinois", "ILLINOIS")
states$ID <- str_replace(states$ID, "indiana", "INDIANA")

states$nudge_y=0
states$nudge_y[states$ID=="ILLINOIS"] <- 1.5
states$nudge_y[states$ID=="INDIANA"] <- 1.6
states$nudge_y[states$ID=="WISCONSIN"] <- -0.8
states$nudge_y[states$ID=="MICHIGAN"] <- -0.5

states$nudge_x=0
states$nudge_x[states$ID=="ILLINOIS"] <- 1
states$nudge_x[states$ID=="WISCONSIN"] <- 1.7
states$nudge_x[states$ID=="MICHIGAN"] <- 0

##color
n=length(unique(loc$pop))
popcol=setNames(pal_igv("default")(17), levels(loc$pop))


p_map <- ggplot(data = world) +
  geom_sf() +
  geom_sf(data = states, fill = "darkgray") + 
  annotate(geom = "text", x = -87, y = 43.5, label = "LAKE MICHIGAN", 
           fontface = "bold", size = 5, angle=70) +
  annotate(geom="text", x=-87.6, y=45, label="Green Bay",
           fontface="bold", size=4, angle=58)+
  annotate(geom="text", x=-87.25, y=45,label="Door County",
           fontface="bold", size=4, angle=58)+
  annotation_scale(location = "br", width_hint = 0.3) +
  geom_text(data = states, aes(X, Y, label = ID), size = 5, fontface="bold", 
            nudge_y=states$nudge_y, nudge_x = states$nudge_x) +
  coord_sf(xlim = c(-88.7, -84.7), ylim = c(41.5, 46.5), expand = TRUE)+
  geom_point(data=loc, aes(x=long, y=lat, color=pop), alpha=0.8, size=2)+
  geom_label_repel(data=loc, aes(x=long, y=lat, label=code, color=pop),
                   size=5, 
                   box.padding = 1.2, max.overlaps = Inf, show.legend = FALSE)+ #show.legend remove legend for text color
  scale_color_manual(values=popcol)+
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.47, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering)+
  labs(x="Longitutde", y="Latitude")+
  theme_bw(base_size=20)+
  theme(legend.position = "none") 

pdf("./figures/fig1A_map.pdf", width =6, height = 12)
p_map
dev.off()



