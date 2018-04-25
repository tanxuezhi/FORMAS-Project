#####################
### Illustrations ###
#####################
source("functions.R")
library(rgeos)

##############
## data set ##
##############

butterflies.data.presence <- as.tbl(fread("../Data/cti_butterflies_data.csv")) %>% dplyr::filter(type == "Presence")
butterflies.data.presence %>% dplyr::filter(Scale == 50000) %>% group_by(country) %>% summarise(n = length(unique(Site)))
butterflies.data.presence %>% dplyr::filter(Scale == 50000) %>% group_by(country) %>% summarise(first = min(Year))



data <- butterflies.data.presence %>% dplyr::filter(Scale == 50000) %>% 
  group_by(Site) %>% summarise(nYear = n(), X = unique(X), Y = unique(Y))

ggplot(data, aes(x = X, y = Y, color = nYear)) + geom_point(alpha = .5) + 
  scale_color_gradientn(name = "No. years \n", 
                       breaks =  c(1,max(data$nYear)),
                       limits = c(1,max(data$nYear)),
                       colours = c("yellow","red")) +
  scale_x_continuous("Longitude") + scale_y_continuous("Latitude")


###############
## landscape ##
###############
ggplot(butterflies.data.presence %>% dplyr:::filter(PLAND < .95 & PLAND > .05) %>% group_by(country, Scale, Habitat) %>%
         distinct(), aes(y = CLUMPY, x = PLAND, color =country)) + 
  geom_point(alpha = .3) + facet_wrap( ~ Scale, labeller = label_wrap) + 
  scale_color_manual("Country", values = c("#00008B", "#BDB76B"), labels = c("Finland", "Netherlands")) + 
  scale_y_continuous("Habitat clumpiness") + scale_x_continuous("% Semi-natural habitat")

SNH_open <- raster("../Landcover/SNH/SNH_open.tif")
SNH_forest <- raster("../Landcover/SNH/SNH_forest.tif")

butterflies.data.presence %>% dplyr::filter(Scale == 20000, CLUMPY < .8, PLAND > .7)

# high clump, high snh
pt1 <- butterflies.data.presence %>% dplyr::filter(Site == "1231_NL") %>% dplyr::select(X, Y) %>% summarise_all(unique)
buf1 <- gBuffer(spgeom = SpatialPoints(pt1), width = 20000)
land1 <- raster::mask(crop(SNH_open, buf1), buf1)

# high clump, low snh
pt2 <- butterflies.data.presence %>% dplyr::filter(Site == "1743_NL") %>% dplyr::select(X, Y) %>% summarise_all(unique)
buf2 <- gBuffer(spgeom = SpatialPoints(pt2), width = 20000)
land2 <- raster::mask(crop(SNH_open, buf2), buf2)

# low clump, high snh
pt3 <- butterflies.data.presence %>% dplyr::filter(Site == "72_FIN") %>% dplyr::select(X, Y) %>% summarise_all(unique)
buf3 <- gBuffer(spgeom = SpatialPoints(pt3), width = 20000)
land3 <- raster::mask(crop(SNH_forest, buf3), buf3)

# low clump, low snh
pt4 <- butterflies.data.presence %>% dplyr::filter(Site == "77_FIN") %>% dplyr::select(X, Y) %>% summarise_all(unique)
buf4 <- gBuffer(spgeom = SpatialPoints(pt4), width = 20000)
land4 <- raster::mask(crop(SNH_open, buf4), buf4)

# plot
par(mfrow=c(2,2))
plot(land2, legend = F, axes = F, box = F)
plot(land1, legend = F, axes = F, box = F)
plot(land4, legend = F, axes = F, box = F)
plot(land3, legend = F, axes = F, box = F)



