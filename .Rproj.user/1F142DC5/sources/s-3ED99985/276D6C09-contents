library(sf)
library(sp)
library(gdistance)
library(raster)
library(tidyverse)
library(rgdal)

rm(list=ls())

####################
###water distance###
####################

## Load Great Lakes shapefile
GL <- st_read("./data/Great_Lakes-shp/Great_Lakes.shp") %>%
  st_transform(crs = 32619)


##lat long
#note that: the lat long info of NMNL and FOXR are located outside of the shape file, 
#therefore, I "corrected" their lat long info so the points are located within the lake;
geo <- read.csv("./data/meta/collections_lat_longs_afterfiltering_corrected.csv") %>% 
  dplyr::select(pop=code, lat, long)


sites <- st_as_sf(geo, coords = c("long", "lat"),
                  crs = 4326, stringsAsFactors = FALSE) %>%
  st_transform(crs = 32619)

rasterGL <- raster(x = extent(GL), nrow = 2000,  ncol = 2000)
rasterGL <- rasterize(x = GL, y = rasterGL, field = 1)
rasterGL[is.na(rasterGL)] <- 0

#create a transition object from the raster
#A transition layer is an object that allows R to estimate least cost paths between two locations. 
#directions: determines the complexity of R's distance tracking techniques.when direction=16, the animal can move to 16 directions;
#It can take a long time to make a transition layer. Save it for a later use; 
rasterGL.TR <- transition(rasterGL, transitionFunction = mean, directions = 16)
#A projection correction is needed for accuracy in the case of grid data for a longlat raster 
rasterGL.TR <- geoCorrection(rasterGL.TR, type = "c")

saveRDS(rasterGL.TR, "./data/meta/GL_transition.rds")

rasterGL.TR <- readRDS("./data/meta/GL_transition.rds")


#the units are undefined since they are a factor of the raster resolution, but it works.
dist <- costDistance(rasterGL.TR,
                      fromCoords = as(as_Spatial(sites), "SpatialPoints"),
                      toCoords = as(as_Spatial(sites), "SpatialPoints"))
colnames(dist) <- geo$pop
rownames(dist) <- geo$pop

write.table(dist, "./data/meta/lwf_waterDist.txt", quote=F, row.names = T)





  
