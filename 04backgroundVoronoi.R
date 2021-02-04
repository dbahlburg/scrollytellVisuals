library(deldir)

myCrs <- crs('+proj=stere +lat_0=-90 +lon_0=-15 +k=1 +x_0=-15 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0')

world <- ne_download(scale = 50, type = 'land',category = 'physical')
antarctica <- world %>% 
  dplyr::mutate(lon = coordinates(world)[,1], lat = coordinates(world)[,2]) %>% 
  dplyr::filter(lat < -54)
antarctica <- spTransform(antarctica, myCrs)


meanRPoly <- rasterToPolygons(meanR, dissolve=TRUE)
hexPoints <- as.data.frame(spsample(meanRPoly, type="random", n=10000)) %>% 
  rowwise() %>% 
  dplyr::mutate(nearestPointInd = which.min((x - abundanceXYZ$x)^2 +
                                            (y - abundanceXYZ$y)^2 ),
                value = abundanceXYZ$meanvalue[nearestPointInd]) %>% 
  ungroup() %>% 
  mutate(weight = value/max(value)) %>% 
  slice(sample(1:10000, size = 2800, prob = value^0.15))
  

try <- deldir(hexPoints$x,hexPoints$y,wl='tr')

tryDf1 <- try[[1]] %>% 
  mutate(distance = sqrt((x2-x1)^2 + (y2-y1)^2)) %>% 
  filter(distance < 250001) %>% 
  rowwise() %>% 
  dplyr::mutate(nearestPointInd1 = which.min((x1 - abundanceXYZ$x)^2 +
                                              (y1 - abundanceXYZ$y)^2 ),
                valueP1 = abundanceXYZ$meanvalue[nearestPointInd1],
                nearestPointInd2 = which.min((x2 - abundanceXYZ$x)^2 +
                                               (y2 - abundanceXYZ$y)^2 ),
                valueP2 = abundanceXYZ$meanvalue[nearestPointInd2],
                weight = valueP1/max(valueP1)) 

lst_result <- data.frame(data = t(as.data.frame(apply(tryDf1, 1, function(x){
  var1 <- as.list(seq(x[['x1']],x[['x2']],length.out = 10))
  var2 <- as.list(seq(x[['y1']],x[['y2']],length.out = 10))
  var3 <- as.list(seq(x[['valueP1']],x[['valueP2']],length.out = 10))
  return(data.frame(xseq = var1, yseq = var2, valseq = var3))
  # return(data.frame(var1 = var1, var2 = var2))
}))))

lst_result2 <- lst_result %>% 
  mutate(pointNo = rep(1:nrow(tryDf1), each = 30),
         pointID = rep(c('x','y','value'), times = nrow(tryDf1), each = 10),
         order = rep(1:10, times = nrow(tryDf1)*3)) %>% 
  pivot_wider(values_from = data, names_from = pointID) 


pp<-lst_result2 %>% 
  ggplot(.,aes(x = x, y = y)) +
  geom_path(aes(alpha = value, size = value, group = pointNo), 
            colour = '#d9685b') +
  scale_size_continuous(range = c(0.6,1)) +
  scale_alpha_continuous(range = c(0.125,1)) +
  geom_point(data = hexPoints, aes(alpha = value, size = value),
             colour = '#ffffff', shape = 16) +
  new_scale('size') +
  new_scale('alpha') +
  scale_size_continuous(range = c(0.05,0.25)) +
  scale_alpha_continuous(range = c(0.5,1)) +
  scale_x_continuous(limits = c(-4e6,4e6)) +
  scale_y_continuous(limits = c(-4.5e6,4.5e6)) +
  # geom_path(colour = '#ffffff') +
  # geom_point(size = 0.1, colour = '#ffffff') +
  geom_polygon(data = antarctica, 
               aes(x = long, y = lat, group = group),
               fill = '#000000', colour = '#ffffff',
               size = 0.2)+
  theme(panel.background = element_rect(fill = 'black', colour = NA),
        plot.background = element_rect(fill = 'black', colour = NA),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none') 
ggsave('~/github/scrollytellJs/visuals/backgroundVoronoi3.png', pp, 
       width = 20, height = 20, bg = 'transparent')





