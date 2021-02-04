library(ncdf4)
library(raster)
library(tidyverse)
library(ggnewscale)
library(deldir)
library(rnaturalearth)

world <- ne_download(scale = 50, type = 'land',category = 'physical')

antarctica <- crop(world, extent(-180,180,-90,-54))
# dplyr::mutate(lon = coordinates(world)[,1], lat = coordinates(world)[,2]) %>% 
# dplyr::filter(lat < -54)
antarctica <- spTransform(antarctica, myCrs)

chl1 <- brick('~/Downloads/dataset-oc-glo-bio-multi-l4-chl_4km_monthly-rep_1611872677115.nc')
chl2 <- brick('~/Downloads/dataset-oc-glo-bio-multi-l4-chl_4km_monthly-rep_1611872916529.nc')
chl3 <- brick('~/Downloads/dataset-oc-glo-bio-multi-l4-chl_4km_monthly-rep_1611872979514.nc')
chl4 <- brick('~/Downloads/dataset-oc-glo-bio-multi-l4-chl_4km_monthly-rep_1611873011531.nc')
chl5 <- brick('~/Downloads/dataset-oc-glo-bio-multi-l4-chl_4km_monthly-rep_1611873041600.nc')
chl6 <- brick('~/Downloads/dataset-oc-glo-bio-multi-l4-chl_4km_monthly-rep_1611873077142.nc')
chl7 <- brick('~/Downloads/dataset-oc-glo-bio-multi-l4-chl_4km_monthly-rep_1611873103366.nc')
chl8 <- brick('~/Downloads/dataset-oc-glo-bio-multi-l4-chl_4km_monthly-rep_1611873132075.nc')
chl9 <- brick('~/Downloads/dataset-oc-glo-bio-multi-l4-chl_4km_monthly-rep_1611873158274.nc')
chl10 <- brick('~/Downloads/dataset-oc-glo-bio-multi-l4-chl_4km_monthly-rep_1611873188363.nc')

allData <- stack(chl1,chl2,chl3,chl4,chl5,chl6,chl7,chl8,chl9,chl10)
allDataAv <- calc(allData, fun = mean, na.rm = T)
myCrs <- crs('+proj=stere +lat_0=-90 +lon_0=-15 +k=1 +x_0=-15 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0')
allDataAvStereo <- projectRaster(allDataAv, crs = myCrs)

southernOceanAvDf <- as.data.frame(allDataAvStereo, xy = T) %>% 
  filter(!is.na(layer)) 

southernOceanAvDfSample <- southernOceanAvDf %>% 
  mutate(layerNorm = layer/max(layer)) %>% 
  slice(sample(1:nrow(southernOceanAvDf), 2200, prob = layerNorm^1)) %>% 
  filter(between(x, -4e6,4e6)) %>% 
  filter(between(y, -4.5e6,4.5e6)) 
southernOceanAvDfSample$quantValue <- ntile(southernOceanAvDfSample$layer, 5)

try <- deldir(southernOceanAvDfSample$x,southernOceanAvDfSample$y,wl='tr')

tryDf1 <- try[[1]]%>% 
  mutate(distance = sqrt((x2-x1)^2 + (y2-y1)^2)) %>% 
  filter(distance < 400001) 

valueP1 <- raster::extract(allDataAvStereo, tibble(x = tryDf1$x1,
                                           y = tryDf1$y1))
valueP2 <- raster::extract(allDataAvStereo, tibble(x = tryDf1$x2,
                                           y = tryDf1$y2))

tryDf1 <- tryDf1 %>% 
  mutate(valueP1 = valueP1,
         valueP2 = valueP2)

lst_result <- data.frame(data = t(as.data.frame(apply(tryDf1, 1, function(x){
  var1 <- as.list(seq(x[['x1']],x[['x2']],length.out = 6))
  var2 <- as.list(seq(x[['y1']],x[['y2']],length.out = 6))
  var3 <- as.list(seq(x[['valueP1']],x[['valueP2']],length.out = 6))
  return(data.frame(xseq = var1, yseq = var2, valseq = var3))
  # return(data.frame(var1 = var1, var2 = var2))
}))))

lst_result2 <- lst_result %>% 
  mutate(pointNo = rep(1:nrow(tryDf1), each = 18),
         pointID = rep(c('x','y','value'), times = nrow(tryDf1), each = 6),
         order = rep(1:6, times = nrow(tryDf1)*3)) %>% 
  pivot_wider(values_from = data, names_from = pointID) 

lst_result2$quantValue <- ntile(lst_result2$value, 5)

pp <- lst_result2 %>% 
  ggplot(.,aes(x = x, y = y)) +
  geom_path(aes(size = log10(value), alpha = log10(value), group = pointNo),
            colour = '#79d491') +#, size = 1.3, alpha = 0.8)+
  scale_size_continuous(range = c(0.02,1.5)) +
  scale_alpha_continuous(range = c(0,0.9)) +
  geom_point(data = southernOceanAvDfSample,
             aes(size = log10(layer), alpha = log10(layer)),
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
ggsave('~/github/scrollytellJs/visuals/primaryProdVoronoi2.png', pp, 
       width = 22, height = 20, bg = 'transparent')




# ggplot(as.data.frame(allDataAvStereo, xy = T), aes(x = x, y = y, fill = log10(layer))) +
#   geom_raster() +scale_fill_viridis_c() +
#   geom_polygon(data = antarctica,
#                aes(x = long, y = lat, group = group),
#                fill = '#000000', colour = '#ffffff',
#                size = 0.15)+
#   theme(panel.background = element_rect(fill = 'black', colour = NA),
#         plot.background = element_rect(fill = 'black', colour = NA),
#         panel.grid = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         legend.position = 'none') 


