library(tidyverse)
library(raster)
library(spdplyr)
library(maptools)
library(rgeos)
library(rnaturalearth)
world <- ne_download(scale = 50, type = 'land',category = 'physical')
antarctica <- world %>% 
  dplyr::mutate(lon = coordinates(world)[,1], lat = coordinates(world)[,2]) %>% 
  dplyr::filter(lat < -60)

setwd("~/Dropbox/NasaTypeAnimation")

allSeaIce <- brick('seaIce/seaIce2017.nc')
seaIce <- subset(allSeaIce,180)
seaIceAgg <- aggregate(seaIce, fact=10)
seaIceAgg <- ceiling(seaIceStereoAgg / 1) * 1
seaIceAggRas <- rasterToPolygons(seaIceAgg, dissolve = T)


seaIceAggRasEdit <- seaIceAggRas %>% 
  dplyr::mutate(lon = coordinates(seaIceAggRas)[,1], lat = coordinates(seaIceAggRas)[,2],
                lat = ifelse(lat < -60, -90, lat))

roughOutline <- tibble(
  lon = c(
    -60.7407622,-78.8868284,-100.6778769,-122.2931912,-143.3813026,-175.1892041,-180,-180,-104.0168279, 28.4868058, 137.6177827, 180, 180, 180, 180, 180, 173.6008765, 163.9354921, 147.5682889, 129.5794438, 100.2318220, 73.4323470, 59.2857389, 45.9299350, 32.2226625, 23.1723480, 8.3228029,-1.3425816,-12.9410429,-21.2005532,-29.6357978,-40.3555878,-49.4937694,-60.7407622
  ),
  lat = c(
    -59.9330004,-64.3208716,-67.1358294,-69.3493386,-69.2872570,-69.5958901,-75.3645057,-84.8183734,-84.5078156,-84.7060489,-84.5580571,-85.6621392,-83.6575540,-81.9724313,-79.8123023,-66.7919095,-66.1605106,-64.6615174,-62.4310742,-61.8561488,-61.6481625,-62.3087937,-62.4717237,-62.9552230,-66.1249624,-66.9300603,-67.0331628,-65.7306265,-64.6991054,-64.2063772,-62.5123179,-59.8889369,-58.9500082,-59.9330004
  ),
  id = 'a'
) 

# make a list
roughOutlineList <- split(roughOutline, roughOutline$id)
# only want lon-lats in the list, not the names
roughOutlineList <- lapply(roughOutlineList, function(x) { x["id"] <- NULL; x })

# create SpatialPolygons Object, convert coords to polygon
ps <- sapply(roughOutlineList, Polygon)

# add id variable 
p1 <- Polygons(ps, ID = 1) 

# create SpatialPolygons object
roughOutlinePoly <- SpatialPolygons(list(p1), proj4string = CRS("+proj=longlat +datum=WGS84") ) 

everything <- union(seaIceAggRas, roughOutlinePoly)
hexPoints <- spsample(everything, type="hexagonal", cellsize=4.5)
hexagons <- fortify(HexPoints2SpatialPolygons(hexPoints))


#assign value to hexagons identifying whether ice or land is represented
hexCenters <- hexagons %>% 
  group_by(group) %>% 
  dplyr::summarise(latHex = mean(lat),
                   longHex = mean(long)) %>% 
  mutate(value =  extract(seaIce, 
                          tibble(long = longHex, lat = latHex))) %>% 
  dplyr::select(-latHex, -longHex) 

hexagonsPlot <- hexagons %>% 
  left_join(hexCenters, by = c('group'))  %>% 
  mutate(value = ifelse(is.na(value), 1, value)) 

ggplot(hexagonsPlot) +
  geom_polygon(aes(x = long, y = lat, group = group, alpha = value),
                fill = '#ff5500', size = 0.25, colour = '#ff5500') +
  scale_alpha_continuous(range = c(0.2,0.4)) +
  # geom_sf(data = world, fill = '#1f1f1f') +
  geom_polygon(data = fortify(antarctica), 
               aes(x = long, y = lat, group = group),
               fill = '#1f1f1f', colour = '#1f1f1f',
               size = 0.1)+
  scale_y_continuous(limits = c(-90,30)) +
  theme(panel.background = element_rect(fill = '#2e2e2e'),
        plot.background = element_rect(fill = '#2e2e2e'),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')

ggsave('background.png', width = 35, height = 12)
