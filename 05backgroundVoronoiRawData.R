krillbaseAgg <- krillbase %>% 
  mutate(Latitude = round(Latitude, 2),
         Longitude = round(Longitude, 2)) %>% 
  dplyr::group_by(Latitude, Longitude) %>% 
  dplyr::summarise(meanBiomass = mean(`Standardisedkrillunder 1m2`, na.rm = T)) %>% 
  filter(!is.na(meanBiomass)) %>% 
  filter(meanBiomass > 0) #%>% 

krillbasePoints <- krillbaseAgg %>% 
  #dplyr::select(-meanBiomass) %>% 
  relocate(Longitude, Latitude)

krillSP <- SpatialPoints(krillbasePoints)
crs(krillSP) <- CRS("+init=epsg:4326")
krillSP <- as.data.frame(spTransform(krillSP, myCrs))
# write.csv(krillSP, '~/Desktop/krillbaseStations.csv')

try <- deldir(krillSP$Longitude,krillSP$Latitude,wl='tr')


tryDf1 <- try[[1]] %>% 
  mutate(euklidDistance = sqrt((x2-x1)^2 + (y2-y1)^2)) %>% 
  filter(euklidDistance < 1000001) %>% 
  rowwise() %>% 
  dplyr::mutate(nearestPointInd1 = which.min((x1 - krillSP$Longitude)^2 +
                                               (y1 - krillSP$Latitude)^2 ),
                valueP1 = krillSP$meanBiomass[nearestPointInd1],
                nearestPointInd2 = which.min((x2 - krillSP$Longitude)^2 +
                                               (y2 - krillSP$Latitude)^2 ),
                valueP2 = krillSP$meanBiomass[nearestPointInd2]) %>% 
  ungroup() %>% 
  mutate(weight = valueP1/max(valueP1)) %>% 
  dplyr::select(x1,y1,x2,y2,ind1,ind2,euklidDistance,valueP1,valueP2)
write.csv(tryDf1, '~/Desktop/krillbaseStationsVertices.csv')


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
            colour = '#95b5e6') +
  scale_size_continuous(range = c(0.9,1.2)) +
  scale_alpha_continuous(range = c(0.25,1)) +
  geom_point(data = krillSP, aes(x = Longitude, 
                                      y = Latitude,
                                      alpha = meanBiomass,
                                      size = meanBiomass),
             colour = '#ffffff', shape = 16) +
  new_scale('size') +
  new_scale('alpha') +
  scale_size_continuous(range = c(0.3,1.3)) +
  scale_alpha_continuous(range = c(0.25,1)) +
  # geom_path(colour = '#ffffff') +
  # geom_point(size = 0.1, colour = '#ffffff') +
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
ggsave('~/github/scrollytellJs/visuals/backgroundVoronoi2.png', pp, 
       width = 20, height = 20, bg = 'transparent')

