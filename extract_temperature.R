library(raster)

######
# NL #
######

# sites data
sites_NL <- read_csv(file = "../Data/Butterflies - Netherlands/Sites_NL_ETRS89_landcover.csv")
sites_NL_WGS84 <- spTransform(SpatialPoints(sites_NL[,c("X", "Y")], proj4string = CRS("+init=epsg:3035")), 
                              CRS("+init=epsg:4326"))

# add mean temperature
temperature.NL <- raster:::getData('worldclim', var='tmean', res = 0.5, 
                                   lon = mean(sites_NL_WGS84@coords[,1]), 
                                   lat = mean(sites_NL_WGS84@coords[,2]),
                                   path = "//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Climate data/")
temperature.NL <- mean(temperature.NL)
temperature.NL.pts <- cbind(sites_NL, temp=extract(temperature.NL, sites_NL_WGS84))

# add NA
sites_NL_NA <- temperature.NL.pts %>% dplyr:::filter(is.na(temp))
sites_NL_NA_WGS84 <- spTransform(SpatialPoints(sites_NL_NA[,c("X", "Y")], proj4string = CRS("+init=epsg:3035")), 
                                  CRS("+init=epsg:4326"))
temperature.NL.large <- aggregate(temperature.NL, 4)
temperature.NL.NA.pts <- cbind(sites_NL_NA[,-5], temp=extract(temperature.NL.large, sites_NL_NA_WGS84))
temperature.NL.pts[temperature.NL.pts$Site %in% (temperature.NL.NA.pts %>% dplyr:::filter(!is.na(temp)))$Site, "temp"] <- (temperature.NL.NA.pts %>% dplyr:::filter(!is.na(temp)))$temp

sites_NL_NA <- temperature.NL.pts %>% dplyr:::filter(is.na(temp))
sites_NL_NA_WGS84 <- spTransform(SpatialPoints(sites_NL_NA[,c("X", "Y")], proj4string = CRS("+init=epsg:3035")), 
                                 CRS("+init=epsg:4326"))
temperature.NL.large <- aggregate(temperature.NL, 8)
temperature.NL.NA.pts <- cbind(sites_NL_NA[,-5], temp=extract(temperature.NL.large, sites_NL_NA_WGS84))
temperature.NL.pts[temperature.NL.pts$Site %in% (temperature.NL.NA.pts %>% dplyr:::filter(!is.na(temp)))$Site, "temp"] <- (temperature.NL.NA.pts %>% dplyr:::filter(!is.na(temp)))$temp

write_csv(temperature.NL.pts %>% dplyr:::select(Site, temp), "../Data/temperature_NL.csv")


#######
# FIN #
#######

# sites data
sites_FIN <- read_csv(file = "../Data/Butterflies - Finland/Sites_FIN_ETRS89_landcover.csv")
sites_FIN_WGS84 <- spTransform(SpatialPoints(sites_FIN[,c("X", "Y")], proj4string = CRS("+init=epsg:3035")), 
                              CRS("+init=epsg:4326"))

# add mean temperature
temperature.FIN_1 <- raster:::getData('worldclim', var='tmean', res = 0.5, 
                                    lon = mean(sites_FIN_WGS84@coords[,1]), 
                                    lat = mean(sites_FIN_WGS84@coords[,2]),
                                    path = "//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Climate data/")
temperature.FIN_1 <- mean(temperature.FIN_1)
temperature.FIN_2 <- raster:::getData('worldclim', var='tmean', res = 0.5, 
                                    lon = max(sites_FIN_WGS84@coords[,1]), 
                                    lat = mean(sites_FIN_WGS84@coords[,2]),
                                    path = "//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Climate data/")
temperature.FIN_2 <- mean(temperature.FIN_2)
temperature.FIN_3 <- raster:::getData('worldclim', var='tmean', res = 0.5, 
                                      lon = mean(sites_FIN_WGS84@coords[,1]), 
                                      lat = min(sites_FIN_WGS84@coords[,2]),
                                      path = "//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Climate data/")
temperature.FIN_3 <- mean(temperature.FIN_3)
temperature.FIN <- merge(temperature.FIN_1, temperature.FIN_2, temperature.FIN_3)
temperature.FIN.pts <- cbind(sites_FIN, temp=extract(temperature.FIN, sites_FIN_WGS84))


# add NA
sites_FIN_NA <- temperature.FIN.pts %>% dplyr:::filter(is.na(temp))
sites_FIN_NA_WGS84 <- spTransform(SpatialPoints(sites_FIN_NA[,c("X", "Y")], proj4string = CRS("+init=epsg:3035")), 
                                  CRS("+init=epsg:4326"))
temperature.FIN.large <- aggregate(temperature.FIN, 4)
temperature.FIN.NA.pts <- cbind(sites_FIN_NA[,-5], temp=extract(temperature.FIN.large, sites_FIN_NA_WGS84))
temperature.FIN.pts[temperature.FIN.pts$Site %in% (temperature.FIN.NA.pts %>% dplyr:::filter(!is.na(temp)))$Site, "temp"] <- (temperature.FIN.NA.pts %>% dplyr:::filter(!is.na(temp)))$temp

write_csv(temperature.FIN.pts %>% dplyr:::select(Site, temp), "../Data/temperature_FIN.csv")


#######
# SWE #
#######

# sites data
sites_SWE <- read_csv(file = "../Data/Birds - Sweden/Sites_SWE_ETRS89_landcover.csv")
sites_SWE_WGS84 <- spTransform(SpatialPoints(sites_SWE[,c("X", "Y")], proj4string = CRS("+init=epsg:3035")), 
                               CRS("+init=epsg:4326"))

# add mean temperature
temperature.SWE_1 <- raster:::getData('worldclim', var='tmean', res = 0.5, 
                                    lon = mean(sites_SWE_WGS84@coords[,1]), 
                                    lat = mean(sites_SWE_WGS84@coords[,2]),
                                    path = "//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Climate data/")
temperature.SWE_1 <- mean(temperature.SWE_1)
temperature.SWE_2 <- raster:::getData('worldclim', var='tmean', res = 0.5, 
                                      lon = mean(sites_SWE_WGS84@coords[,1]), 
                                      lat = min(sites_SWE_WGS84@coords[,2]),
                                      path = "//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Climate data/")
temperature.SWE_2 <- mean(temperature.SWE_2)
temperature.SWE <- merge(temperature.SWE_1, temperature.SWE_2)
plot(temperature.SWE);points(sites_SWE_WGS84)
temperature.SWE.pts <- cbind(sites_SWE, temp=extract(temperature.SWE, sites_SWE_WGS84))

# add NA
sites_SWE_NA <- temperature.SWE.pts %>% dplyr:::filter(is.na(temp))
sites_SWE_NA_WGS84 <- spTransform(SpatialPoints(sites_SWE_NA[,c("X", "Y")], proj4string = CRS("+init=epsg:3035")), 
                               CRS("+init=epsg:4326"))
temperature.SWE.large <- aggregate(temperature.SWE, 4)
temperature.SWE.NA.pts <- cbind(sites_SWE_NA[,-5], temp=extract(temperature.SWE.large, sites_SWE_NA_WGS84))
temperature.SWE.pts[temperature.SWE.pts$Site %in% (temperature.SWE.NA.pts %>% dplyr:::filter(!is.na(temp)))$Site, "temp"] <- (temperature.SWE.NA.pts %>% dplyr:::filter(!is.na(temp)))$temp

sites_SWE_NA <- temperature.SWE.pts %>% dplyr:::filter(is.na(temp))
sites_SWE_NA_WGS84 <- spTransform(SpatialPoints(sites_SWE_NA[,c("X", "Y")], proj4string = CRS("+init=epsg:3035")), 
                                  CRS("+init=epsg:4326"))
temperature.SWE.large <- aggregate(temperature.SWE, 90)
temperature.SWE.NA.pts <- cbind(sites_SWE_NA[,-5], temp=extract(temperature.SWE.large, sites_SWE_NA_WGS84))
temperature.SWE.pts[temperature.SWE.pts$Site %in% (temperature.SWE.NA.pts %>% dplyr:::filter(!is.na(temp)))$Site, "temp"] <- (temperature.SWE.NA.pts %>% dplyr:::filter(!is.na(temp)))$temp

write_csv(temperature.SWE.pts %>% dplyr:::select(Site, temp), "../Data/temperature_SWE.csv")

