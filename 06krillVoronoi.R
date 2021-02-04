library(png)
library(tidyverse)
library(deldir)
library(raster)
library(ggnewscale)

krill <- readPNG('~/github/scrollytellVisuals/krillRaster.png')
krillRas <- raster('~/github/scrollytellVisuals/krillRaster2.png')
krillRas[krillRas == 255] <- NA

krillBW <- as_tibble(256 * (krill[,,1] +krill[,,2] +krill[,,3])/3) %>% 
  gather(variable, value) %>% 
  mutate(x = as.numeric(str_replace(variable, pattern = 'V', '')),
         y = rep(dim(krill)[1]:1, times = dim(krill)[2])) %>% 
  filter(value < 256) %>% 
  mutate(value = 1-value/max(value)) %>% 
  slice(sample(1:nrow(.), 400, prob = value^0.8))

# krillBW %>%
#   ggplot(.,aes(x = x, y = y, fill = value)) +
#   geom_raster()
#   

try <- deldir(krillBW$x,krillBW$y,wl='tr')


tryDf1 <- try[[1]] %>% 
  mutate(euklidDistance = sqrt((x2-x1)^2 + (y2-y1)^2),
         centerPointX = x1 + (x2-x1)/2,
         centerPointY = y1 + (y2-y1)/2)

krillMaskValues <- extract(krillRas, tibble(x = tryDf1$centerPointX,
                                               y = tryDf1$centerPointY))
tryDf1 <- tryDf1 %>% 
  mutate(krillMask = krillMaskValues) %>% 
  filter(!is.na(krillMask))
  
krillDf <- tibble(x = c(tryDf1$x1,tryDf1$x2),
                  y = c(tryDf1$y1,tryDf1$y2),
                  ID = rep(1:nrow(tryDf1), times = 2)) %>% 
  left_join(.,krillBW)

pp <- krillDf %>% 
  ggplot(.,aes(x = x, y = y, group = ID)) +
  geom_path(aes(size = value, alpha = value),colour = '#eb6657') +#, size = 1.3, alpha = 0.8)+
  scale_size_continuous(range = c(2.2,2.5)) +
  scale_alpha_continuous(range = c(0.7,0.9)) +
  geom_point(aes(size = value, alpha = value),colour = '#ffffff', shape = 16, size = 3)+#, size = 1.7) +
  #new_scale('size') +
  #scale_size_continuous(range = c(8,10)) +
  theme(panel.background = element_rect(fill = NA, colour = NA),
        plot.background = element_rect(fill = NA, colour = NA),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none') 
ggsave('~/github/scrollytellJs/visuals/krillVoronoi2.png', pp, 
       width = 30, height = 19, bg = 'transparent')

