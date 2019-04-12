library(raster)
library(tidyverse)

### load data
butterflies.data <- bind_rows("FIN" = read_csv(file = "../Data/Butterflies - Finland/Sites_FIN_ETRS89_landcover.csv") %>% 
                                mutate(Site = paste0(Site, "_FIN")),
                              "NL" = read_csv(file = "../Data/Butterflies - Netherlands/Sites_NL_ETRS89_landcover.csv") %>% 
                                mutate(Site = paste0(Site, "_NL")), .id = "country")

birds.data <- read_csv("../Data/Birds - Sweden/Sites_SWE_ETRS89_landcover.csv")

source("fragstat.R")

### load rasters
# #separate forest and open habitat
CLC_FIN <- raster("../Landcover/Corine_land-cover_2012_raster/CLC_FIN.tif")
CLC_NL <- raster("../Landcover/Corine_land-cover_2012_raster/CLC_NL.tif")
CLC_SWE <- raster("../Landcover/Corine_land-cover_2012_raster/CLC_SWE.tif")

CLC <- merge(CLC_FIN, CLC_NL, CLC_SWE)

LC_reclass_open <- as.tbl(rio::import(file = "../Landcover/clc_legend.xls", which = 5L))
LC_reclass_forest<- as.tbl(rio::import(file = "../Landcover/clc_legend.xls", which = 6L))

SNH_open <- reclassify(CLC, cbind.data.frame(is = LC_reclass_open[,1], becomes = LC_reclass_open[,4]))
SNH_forest <- reclassify(CLC, cbind.data.frame(is = LC_reclass_forest[,1], becomes = LC_reclass_forest[,4]))

writeRaster(SNH_open, "../Landcover/SNH/SNH_open.tif", overwrite = T)
writeRaster(SNH_forest, "../Landcover/SNH/SNH_forest.tif", overwrite = T)

SNH_open <- raster("../Landcover/SNH/SNH_open.tif")
SNH_forest <- raster("../Landcover/SNH/SNH_forest.tif")


#one SNH category
CLC_SNH <- raster("../Landcover/SNH/CLC_SNH_cropped.tif")
#crop SNH raster
sites_FIN <- SpatialPointsDataFrame(butterflies.data %>% dplyr::filter(country == "FIN") %>% dplyr::select(X,Y), 
                                           butterflies.data %>% dplyr::filter(country == "FIN"), proj4string = CRS("+init=epsg:3035"))
sites_NL <- SpatialPointsDataFrame(butterflies.data %>% dplyr::filter(country == "NL") %>% dplyr::select(X,Y), 
                                    butterflies.data %>% dplyr::filter(country == "NL"), proj4string = CRS("+init=epsg:3035"))
sites_SWE <- SpatialPointsDataFrame(birds.data %>% dplyr::select(X,Y), 
                                    birds.data, proj4string = CRS("+init=epsg:3035"))

CLC_SNH_FIN <- crop(CLC_SNH, extent(sites_FIN) + 100000)
CLC_SNH_NL <- crop(CLC_SNH, extent(sites_NL) + 100000)
CLC_SNH_SWE <- crop(CLC_SNH, extent(sites_SWE) + 100000)

CLC_SNH_cropped <- merge(CLC_SNH_FIN, CLC_SNH_NL, CLC_SNH_SWE)
writeRaster(CLC_SNH_cropped, "../Landcover/SNH/CLC_SNH_cropped.tif", overwrite = T)

CLC_SNH_cropped <- raster("../Landcover/SNH/CLC_SNH_cropped.tif")
# 
# #write Fragstat files
# butterflies.data <- butterflies.data %>% group_by(coords = paste(X,Y)) %>% summarise(X = unique(X), Y = unique(Y),
#                                                                                      country = unique(country),
#                                                                                      Landcover = unique(Landcover)) %>% 
#   mutate(LID = 1:nrow(.))
# 
# butterflies.data %>% 
#   mutate(row = rowFromY(CLC_SNH_cropped, Y), col = colFromX(CLC_SNH_cropped, X)) %>% dplyr::select(LID, row, col) %>%
#   with(data.frame(FPT_TABLE = paste0("[",LID,":",row,":",col,"]"))) %>%
#   with(write.table(.,"../Connectivity/Fragmentation/sites.fpt", row.names = F, quote = F))
# 
# 

#### use function ####
# land-cover type by habitat
sites_FIN_forest <- SpatialPointsDataFrame(butterflies.data %>% dplyr::filter(country == "FIN", Landcover %in% c(23,24,25)) %>% dplyr::select(X,Y), 
                                           butterflies.data %>% dplyr::filter(country == "FIN", Landcover %in% c(23,24,25)), proj4string = CRS("+init=epsg:3035"))
sites_FIN_open <- SpatialPointsDataFrame(butterflies.data %>% dplyr::filter(country == "FIN", !Landcover %in% c(23,24,25)) %>% dplyr::select(X,Y), 
                                         butterflies.data %>% dplyr::filter(country == "FIN", !Landcover %in% c(23,24,25)), proj4string = CRS("+init=epsg:3035"))

sites_NL_forest <- SpatialPointsDataFrame(butterflies.data %>% dplyr::filter(country == "NL", Landcover %in% c(311,312,313)) %>% dplyr::select(X,Y), 
                                          butterflies.data %>% dplyr::filter(country == "NL", Landcover %in% c(311,312,313)), proj4string = CRS("+init=epsg:3035"))
sites_NL_open <- SpatialPointsDataFrame(butterflies.data %>% dplyr::filter(country == "NL", !Landcover %in% c(311,312,313)) %>% dplyr::select(X,Y), 
                                        butterflies.data %>% dplyr::filter(country == "NL", !Landcover %in% c(311,312,313)), proj4string = CRS("+init=epsg:3035"))

sites_SWE_forest <- SpatialPointsDataFrame(birds.data %>% dplyr::filter(Landcover %in% c(23,24,25)) %>% dplyr::select(X,Y), 
                                           birds.data %>% dplyr::filter(Landcover %in% c(23,24,25)), proj4string = CRS("+init=epsg:3035"))
sites_SWE_open <- SpatialPointsDataFrame(birds.data %>% dplyr::filter(!Landcover %in% c(23,24,25)) %>% dplyr::select(X,Y), 
                                         birds.data %>% dplyr::filter(!Landcover %in% c(23,24,25)), proj4string = CRS("+init=epsg:3035"))

SNH_open <- "../Landcover/SNH/SNH_open.tif"
SNH_forest <- "../Landcover/SNH/SNH_forest.tif"


library(foreach)
scale = c(1000,3000,5000,10000,20000,30000,50000)
sites = c(sites_SWE_open)

res_frag_open <- foreach (i = scale, .combine = rbind) %:% foreach(j = sites, .combine = rbind) %do% {
  cbind.data.frame(Scale = i, fragstat(points = j, LC_SNH = SNH_open, width = i, cores = 2))
}

sites = c(sites_SWE_forest)

res_frag_forest <- foreach (i = scale, .combine = rbind) %:% foreach(j = sites, .combine = rbind) %do% {
  cbind.data.frame(Scale = i, fragstat(points = j, LC_SNH = SNH_forest, width = i, cores = 2))
}


res_frag <- bind_rows("Open" = res_frag_open, "Forest" = res_frag_forest, .id = "Habitat")
res_frag <- res_frag %>% mutate(country = ifelse(grepl("FIN", .$Site), "FIN", "NL"))

table(res_frag[,-c(3:5)])

ggplot(res_frag %>% dplyr::filter(PLAND > .01), aes(x = PLAND, y = CLUMPY, color = Habitat)) + geom_point() + facet_wrap( ~ Scale)


## write results ##
write_csv(res_frag, "../Connectivity/Fragmentation/Frag_indices2.csv")

# 3 land-cover type for all sites
sites_FIN <- SpatialPointsDataFrame(butterflies.data %>% dplyr::filter(country == "FIN") %>% dplyr::select(X,Y), 
                                    butterflies.data %>% dplyr::filter(country == "FIN"), proj4string = CRS("+init=epsg:3035"))
sites_NL <- SpatialPointsDataFrame(butterflies.data %>% dplyr::filter(country == "NL") %>% dplyr::select(X,Y), 
                                   butterflies.data %>% dplyr::filter(country == "NL"), proj4string = CRS("+init=epsg:3035"))

# generalist (SNH as open SNH + forest)
SNH <- "../Landcover/SNH/CLC_SNH_cropped.tif"

scale = c(1000,3000,5000,10000,20000,30000,50000)
sites = c(sites_FIN, sites_NL)

res_frag_gen <- foreach (i = scale, .combine = rbind) %:% foreach(j = sites, .combine = rbind) %do% {
  cbind.data.frame(Scale = i, fragstat(points = j, LC_SNH = SNH, width = i, cores = 4))
}

# forest
res_frag_forest <- foreach (i = scale, .combine = rbind) %:% foreach(j = sites, .combine = rbind) %do% {
  cbind.data.frame(Scale = i, fragstat(points = j, LC_SNH = SNH_forest, width = i, cores = 4))
}

# open
res_frag_open <- foreach (i = scale, .combine = rbind) %:% foreach(j = sites, .combine = rbind) %do% {
  cbind.data.frame(Scale = i, fragstat(points = j, LC_SNH = SNH_open, width = i, cores = 4))
}


# write
res_frag_Allhab <- bind_rows("Generalist" = res_frag_gen,
                      "Open" = res_frag_open,
                      "Forest" = res_frag_forest, .id = "Habitat")
res_frag_Allhab <- res_frag_Allhab %>% mutate(country = ifelse(grepl("FIN", .$Site), "FIN", "NL"))

table(res_frag_Allhab[,-c(3:5)])

write_csv(res_frag_Allhab, "../Connectivity/Fragmentation/Frag_indices_Allhab.csv")


# 
# ## load Fragstat results
# res_frag <- bind_rows("3000" = read_csv("../Connectivity/Fragmentation/frag_3000.class", col_types = "ccdddddddddd"),
#                       "30000" = read_csv("../Connectivity/Fragmentation/frag_30000.class", col_types = "ccdddddddddd"), .id = "Scale") %>% 
#   mutate(LID = as.numeric(gsub("point_", "", .$LID))) %>% dplyr::select(-TYPE, -PAFRAC)
# 
# left_join(res_frag, butterflies.data) %>% summarise(n = length(unique(LID)))
# left_join(res_frag, butterflies.data) %>% summarise(n = length(unique(paste(X, Y))))
# 
# res_frag <- right_join(butterflies.data  %>% dplyr::select(-Landcover) ,
#                        res_frag, by = "LID") %>% dplyr::select(-LID)
# 
# res_frag %>% group_by(country, Scale) %>% summarise(n = length(unique(coords)))
# 
# ## write Fragstat results
# write_csv(res_frag, "../Connectivity/Fragmentation/Frag_indices.csv")
# 
# ## descriptors ##
# 
# ggplot(res_frag %>% dplyr::filter(PLAND < 95, PLAND > 5), aes(x = PLAND, y = CLUMPY)) + 
#   geom_point() +
#   facet_wrap(Scale ~ country) +
#   scale_y_continuous("Clumping") +
#   scale_x_continuous("% SNH")
