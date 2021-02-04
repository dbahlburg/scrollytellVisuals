library(raster)
library(tidyverse)
library(deldir)
library(rnaturalearth)

world <- ne_download(scale = 50, type = 'land',category = 'physical')

antarctica <- crop(world, extent(-180,180,-90,-54))
antarctica <- spTransform(antarctica, myCrs)


currentsU <- brick('~/Downloads/global-reanalysis-phy-001-031-grepv2-monthly_1611959043931.nc',
                  varname = 'uo_cglo')
currentsV <- brick('~/Downloads/global-reanalysis-phy-001-031-grepv2-monthly_1611959043931.nc',
                  varname = 'vo_cglo')

currentsUAv <- calc(currentsU, fun = mean, na.rm = T)
currentsVAv <- calc(currentsV, fun = mean, na.rm = T)
myCrs <- crs('+proj=stere +lat_0=-90 +lon_0=-15 +k=1 +x_0=-15 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0')
currentsUAvStereo <- projectRaster(currentsUAv, crs = myCrs)
currentsVAvStereo <- projectRaster(currentsVAv, crs = myCrs)

currentsStereo <- sqrt(currentsUAvStereo^2 + currentsVAvStereo^2)
currentsStereoDf <- as.data.frame(currentsStereo, xy=T) %>% 
  filter(!is.na(layer))


currentsStereoDfSample <- currentsStereoDf %>% 
  mutate(layerNorm = layer/max(layer)) %>% 
  slice(sample(1:nrow(currentsStereoDf), 3000, prob = layerNorm^1)) %>% 
  filter(between(x, -4e6,4e6)) %>% 
  filter(between(y, -4.5e6,4.5e6)) 

try <- deldir(currentsStereoDfSample$x,currentsStereoDfSample$y,wl='tr')

tryDf1 <- try[[1]]%>% 
  mutate(distance = sqrt((x2-x1)^2 + (y2-y1)^2)) %>% 
  filter(distance < 500001) 

valueP1 <- raster::extract(currentsStereo, tibble(x = tryDf1$x1,
                                                   y = tryDf1$y1))
valueP2 <- raster::extract(currentsStereo, tibble(x = tryDf1$x2,
                                                   y = tryDf1$y2))

tryDf1 <- tryDf1 %>% 
  mutate(valueP1 = valueP1,
         valueP2 = valueP2)

lst_result <- data.frame(data = t(as.data.frame(apply(tryDf1, 1, function(x){
  var1 <- as.list(seq(x[['x1']],x[['x2']],length.out = 6))
  var2 <- as.list(seq(x[['y1']],x[['y2']],length.out = 6))
  var3 <- as.list(seq(x[['valueP1']],x[['valueP2']],length.out = 6))
  return(data.frame(xseq = var1, yseq = var2, valseq = var3))
}))))

lst_result2 <- lst_result %>% 
  mutate(pointNo = rep(1:nrow(tryDf1), each = 18),
         pointID = rep(c('x','y','value'), times = nrow(tryDf1), each = 6),
         order = rep(1:6, times = nrow(tryDf1)*3)) %>% 
  pivot_wider(values_from = data, names_from = pointID) 


pp <- lst_result2 %>% 
  ggplot(.,aes(x = x, y = y)) +
  geom_path(aes(size = value, alpha = value, group = pointNo),
            colour = '#f75252') +#, size = 1.3, alpha = 0.8)+
  scale_size_continuous(range = c(0.1,1.8)) +
  scale_alpha_continuous(range = c(0.1,1)) +
  geom_point(data = currentsStereoDfSample,
             aes(size = layer, alpha = layer),
             colour = '#ffffff', shape = 16)+#, size = 1.7) +
  new_scale('size') +
  scale_size_continuous(range = c(0.2,1.7)) +
  geom_polygon(data = antarctica,
               aes(x = long, y = lat, group = group),
               fill = '#000000', colour = '#ffffff',
               size = 0.2)+
  scale_x_continuous(limits = c(-4e6,4e6)) +
  scale_y_continuous(limits = c(-4.5e6,4.5e6)) +
  theme(panel.background = element_rect(fill = 'black', colour = NA),
        plot.background = element_rect(fill = 'black', colour = NA),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none') 
ggsave('~/github/scrollytellJs/visuals/advectionVoronoi.png', pp, 
       width = 22, height = 20, bg = 'transparent')








