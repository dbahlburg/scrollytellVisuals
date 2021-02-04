library(tidyverse)
library(raster)
library(spdplyr)
library(maptools)
library(mapproj)
library(rgeos)
library(rnaturalearth)
library(plyr)
library(tidyverse)
library(rnaturalearth)
library(raster)
library(gstat)

krillbase <- read_csv("~/Downloads/krillbase.csv", 
                      col_types = cols(Date = col_character()))


krillbase %>% 
  ggplot(.,aes(x = Longitude, y = Latitude)) +
  geom_point() +
  coord_polar() +
  scale_y_continuous(limits = c(-90, -40))

krillbaseAgg <- krillbase %>% 
  mutate(Latitude = round(Latitude, 2),
         Longitude = round(Longitude, 2)) %>% 
  dplyr::group_by(Latitude, Longitude) %>% 
  dplyr::summarise(meanBiomass = mean(`Standardisedkrillunder 1m2`, na.rm = T)) %>% 
  filter(!is.na(meanBiomass)) %>% 
  filter(meanBiomass > 0)

# points to raster from scratch
df = as.data.frame(krillbaseAgg)
df$cutX = cut(df$Longitude, breaks = seq(-180, 180, 10))
df$cutY = cut(df$Latitude, breaks = seq(-78, -34, 2))

mVal = ddply(df, .(cutX, cutY), summarize, meanvalue = mean(meanBiomass))
mVal$cutX = as.character(mVal$cutX)
mVal$cutX = unlist(strsplit(mVal$cutX, ","))[c(T,F)]
mVal$cutX = as.numeric(sub("\\(", "", mVal$cutX))
mVal$cutY = as.character(mVal$cutY)
mVal$cutY = unlist(strsplit(mVal$cutY, ","))[c(T,F)]
mVal$cutY = as.numeric(sub("\\(", "", mVal$cutY))

coordinates(mVal) = ~cutX +cutY
gridded(mVal) = T
meanr = raster(mVal)
crs(meanr) <- CRS("+init=epsg:4326")
myCrs <- crs('+proj=stere +lat_0=-90 +lon_0=-15 +k=1 +x_0=-15 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0')
meanR = projectRaster(meanr, crs = myCrs)

world <- ne_download(scale = 50, type = 'land',category = 'physical')
antarctica <- crop(world, extent(-180,180,-90,-54))
  # dplyr::mutate(lon = coordinates(world)[,1], lat = coordinates(world)[,2]) %>% 
  # dplyr::filter(lat < -54)
antarctica <- spTransform(antarctica, myCrs)


abundanceXYZ <- as.data.frame(rasterToPoints(meanR))

meanRPoly <- rasterToPolygons(meanR, dissolve=TRUE)
hexPoints <- spsample(meanRPoly, type="hexagonal", cellsize=250000)
hexagons <- fortify(HexPoints2SpatialPolygons(hexPoints))


#assign value to hexagons identifying whether ice or land is represented
hexCenters <- hexagons %>% 
  dplyr::group_by(group) %>% 
  dplyr::summarise(latHex = mean(lat),
                   longHex = mean(long)) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(nearestPointInd = which.min((longHex - abundanceXYZ$x)^2 +
                                       (latHex - abundanceXYZ$y)^2 ),
         value = abundanceXYZ$meanvalue[nearestPointInd]) %>% 
  dplyr::select(-latHex, -longHex, -nearestPointInd) 

hexagonsPlot <- hexagons %>% 
  left_join(hexCenters, by = c('group'))  #%>% 
  filter(value > 5)

pp <- ggplot(hexagonsPlot) +
  geom_polygon(aes(x = long, y = lat, group = group, alpha = value),
               fill = '#ff5500', size = 0.7, colour = '#ff550050') +
  scale_alpha_continuous(range = c(0,0.8)) +
  # geom_sf(data = world, fill = '#1f1f1f') ++
  geom_polygon(data = antarctica, 
               aes(x = long, y = lat, group = group),
               fill = '#1f1f1f', colour = '#ffffff69',
               size = 0.7)+
  theme(panel.background = element_rect(fill = 'transparent', colour = NA),
        plot.background = element_rect(fill = 'transparent', colour = NA),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none') 
ggsave('~/github/scrollytellJs/visuals/backgroundStereo.png', pp, 
       width = 20, height = 20, bg = 'transparent')


#Create legend
hexagonsLegend <- hexagons %>% 
  slice(49:83) %>% 
  mutate(value = rep(c(0:4), each = 7))

legendPlot <- hexagonsLegend %>% 
  ggplot(.,aes(x = long, y = lat, group = group, alpha = value)) +
  geom_polygon(fill = '#ff5500', size = 1, colour = '#ff550050')+
  scale_alpha_continuous(range = c(0,0.8)) +
  scale_x_continuous(limits = c(-2200000, -830000)) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none') 
ggsave('~/github/scrollytellJs/visuals/backgroundStereoLegend.png', legendPlot, 
       width = 10, height = 2.2,bg = "transparent")





