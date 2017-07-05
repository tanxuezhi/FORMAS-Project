library(sp)
library(raster)

CLC <- raster("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Landcover/Corine_land-cover_2012_raster/g100_clc12_V18_5.tif")

### NL ###
Sites <- rio::import(file = "//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Data/Butterflies - Netherlands/Sites.xlsx", which = 1L)
Sites[,-1] <- Sites[,-1]*1000

Sites_RD_New <- SpatialPoints(Sites[,-1], proj4string = CRS("+init=epsg:28992"))
Sites_WGS84 <- spTransform(Sites_RD_New, CRS("+init=epsg:4326"))
Sites_ETRS89 <- spTransform(Sites_RD_New, CRS("+init=epsg:3035"))

write.csv(cbind.data.frame(Site = Sites[,1], 
                           X = Sites_ETRS89@coords[,1], 
                           Y = Sites_ETRS89@coords[,2]),
          "Sites_ETRS89.csv", row.names = F)

### Fin ###
Sites1 <- read.table("../Data/Butterflies - Finland/FINLAND_Sites_WGS84_dgs_1999-2015.txt", sep = "\t", h=T, dec=",")
Sites2 <- read.table("../Data/Butterflies - Finland/FINLAND_Sites_WGS84_dgs_2016.txt", sep = "\t", h=T, dec=",")

Sites <- rbind(Sites1, Sites2)
Sites <- Sites[!duplicated(Sites$Site),]

Sites_WGS84 <- SpatialPoints(Sites[c(5,4)], proj4string = CRS("+init=epsg:4326"))
Sites_ETRS89 <- spTransform(Sites_WGS84, CRS("+init=epsg:3035"))

#add landcover
Landcover <- extract(CLC, Sites_ETRS89)

write.csv(cbind.data.frame(Site = Sites[,1], 
                           X = Sites_ETRS89@coords[,1], 
                           Y = Sites_ETRS89@coords[,2],
                           Landcover = Landcover),
          "../Data/Butterflies - Finland/Sites_FIN_ETRS89_landcover.csv", row.names = F)

