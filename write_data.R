library(rio)
library(sp)
library(raster)
library(tidyverse)

##############################################
##############################################
######                                  ######
###### Community-temperature index data ######
######                                  ######
##############################################
##############################################

lc_class <- as.tbl(rio::import(file = "../Landcover/Corine_land-cover_2012_raster/clc_legend.xls", which = 1L))
lc_class$CLC_CODE <- as.integer(lc_class$CLC_CODE)

###################
### butterflies ###
###################

# sites data
sites_NL <- read_csv(file = "../Data/Butterflies - Netherlands/Sites_NL_ETRS89_landcover.csv")
sites_FIN <- read_csv(file = "../Data/Butterflies - Finland/Sites_FIN_ETRS89_landcover.csv")

# fragmentation data
frag_data <- read_csv("../Connectivity - Fragmentation/Fragmentation/Frag_indices2.csv")

# temperature data
temp_data_NL <- read_csv("../Data/temperature_NL.csv")
temp_data_FIN <- read_csv("../Data/temperature_FIN.csv")

#dbMEM data
# dbMEM_data_NL <- read_csv("../Data/Butterflies - Netherlands/dbMEM.NL.csv")
# dbMEM_data_FIN <- read_csv("../Data/Butterflies - Finland/dbMEM.FIN.csv")

sites_NL <- left_join(sites_NL, lc_class, by = c("Landcover" = "CLC_CODE"))
sites_NL <- left_join(sites_NL, temp_data_NL, by = c("Site"))
sites_NL <- sites_NL %>% dplyr:::select(-GRID_CODE,-RGB)

sites_FIN <- left_join(sites_FIN, lc_class, by = c("Landcover" = "GRID_CODE"))
sites_FIN <- left_join(sites_FIN, temp_data_FIN, by = c("Site"))
sites_FIN <- sites_FIN %>% dplyr:::select(-CLC_CODE,-RGB)

# cti data
# Abundance-weighted #
cti_AB_NL <- read.csv2("../Data/Butterflies - Netherlands/cti_site_year_abundance_2017-05-19.csv", dec = ".")[,-3]
cti_AB_NL <- merge(cti_AB_NL, sites_NL, by.x = "Site", by.y = "Site")

cti_AB_FIN <- read.csv("../Data/Butterflies - Finland/CTI_Abundance_FINLAND_1999-2016.csv", dec = ".")
cti_AB_FIN <- merge(cti_AB_FIN, sites_FIN, by.x = "Site", by.y = "Site")

cti_AB_butterflies <- as.tbl(bind_rows(cbind.data.frame(cti_AB_NL, country = "NL"), cbind.data.frame(cti_AB_FIN, country = "FIN")))
cti_AB_butterflies$Site <- paste0(cti_AB_butterflies$Site, "_", cti_AB_butterflies$country)

# P/A #
cti_PA_NL <- read.csv2("../Data/Butterflies - Netherlands/cti_site_year_presence_2017-05-19.csv", dec = ".")[,-3]
cti_PA_NL <- merge(cti_PA_NL, sites_NL, by.x = "Site", by.y = "Site")

cti_PA_FIN <- read.csv("../Data/Butterflies - Finland/CTI_presence_FINLAND_1999-2016.csv", dec = ".")
cti_PA_FIN <- merge(cti_PA_FIN, sites_FIN, by.x = "Site", by.y = "Site")

cti_PA_butterflies <- as.tbl(bind_rows(cbind.data.frame(cti_PA_NL, country = "NL"), cbind.data.frame(cti_PA_FIN, country = "FIN")))
cti_PA_butterflies$Site <- paste0(cti_PA_butterflies$Site, "_", cti_PA_butterflies$country)

# add grid cell (50x50km)
but.pts <- SpatialPoints(cti_AB_butterflies[,c("X", "Y")])
but.r <- raster(ext = extent(but.pts)*1.1, resolution = 50000)
values(but.r) <- c(1:ncell(but.r))
gridCell50 <- raster::extract(but.r, but.pts)

cti_AB_butterflies <- bind_cols(cti_AB_butterflies, gridCell50 = gridCell50)
cti_PA_butterflies <- bind_cols(cti_PA_butterflies, gridCell50 = gridCell50)

cti_butterflies <- as.tbl(bind_rows(cbind.data.frame(type = "Presence", cti_PA_butterflies), cbind.data.frame(type = "Abundance", cti_AB_butterflies)))

# add fragmentation
cti_butterflies <- left_join(cti_butterflies %>% ungroup(), frag_data) 

cti_butterflies <- cti_butterflies %>% group_by(coords = paste(X, Y), Year, type, Scale) %>% summarise_all(first)

cti_butterflies %>% group_by(Site, Year, type, Scale) %>% summarise(n = n()) %>% ungroup() %>% summarise(max(n))
cti_butterflies %>% group_by(coords, Year, type, Scale) %>% summarise(n = n()) %>% ungroup() %>% summarise(max(n))

cti_butterflies <- cti_butterflies %>% ungroup()  %>% dplyr::select(-coords)

# merge and write
write_csv(cti_butterflies, "../Data/cti_butterflies_data.csv")



#######################
###  alternatively  ###
### same definition ###
###  for all sites  ###
#######################

# sites data
sites_NL <- read_csv(file = "../Data/Butterflies - Netherlands/Sites_NL_ETRS89_landcover.csv")
sites_FIN <- read_csv(file = "../Data/Butterflies - Finland/Sites_FIN_ETRS89_landcover.csv")

# fragmentation data
frag_data <- read_csv("../Connectivity/Fragmentation/Frag_indices_Allhab.csv") %>% filter(Habitat == "Generalist")

# temperature data
temp_data_NL <- read_csv("../Data/temperature_NL.csv")
temp_data_FIN <- read_csv("../Data/temperature_FIN.csv")

#dbMEM data
# dbMEM_data_NL <- read_csv("../Data/Butterflies - Netherlands/dbMEM.NL.csv")
# dbMEM_data_FIN <- read_csv("../Data/Butterflies - Finland/dbMEM.FIN.csv")

sites_NL <- left_join(sites_NL, lc_class, by = c("Landcover" = "CLC_CODE"))
sites_NL <- left_join(sites_NL, temp_data_NL, by = c("Site"))
sites_NL <- sites_NL %>% dplyr:::select(-GRID_CODE,-RGB)

sites_FIN <- left_join(sites_FIN, lc_class, by = c("Landcover" = "GRID_CODE"))
sites_FIN <- left_join(sites_FIN, temp_data_FIN, by = c("Site"))
sites_FIN <- sites_FIN %>% dplyr:::select(-CLC_CODE,-RGB)

# cti data
# Abundance-weighted #
cti_AB_NL <- read.csv2("../Data/Butterflies - Netherlands/cti_site_year_abundance_2017-05-19.csv", dec = ".")[,-3]
cti_AB_NL <- merge(cti_AB_NL, sites_NL, by.x = "Site", by.y = "Site")

cti_AB_FIN <- read.csv("../Data/Butterflies - Finland/CTI_Abundance_FINLAND_1999-2016.csv", dec = ".")
cti_AB_FIN <- merge(cti_AB_FIN, sites_FIN, by.x = "Site", by.y = "Site")

cti_AB_butterflies <- as.tbl(bind_rows(cbind.data.frame(cti_AB_NL, country = "NL"), cbind.data.frame(cti_AB_FIN, country = "FIN")))
cti_AB_butterflies$Site <- paste0(cti_AB_butterflies$Site, "_", cti_AB_butterflies$country)

# P/A #
cti_PA_NL <- read.csv2("../Data/Butterflies - Netherlands/cti_site_year_presence_2017-05-19.csv", dec = ".")[,-3]
cti_PA_NL <- merge(cti_PA_NL, sites_NL, by.x = "Site", by.y = "Site")

cti_PA_FIN <- read.csv("../Data/Butterflies - Finland/CTI_presence_FINLAND_1999-2016.csv", dec = ".")
cti_PA_FIN <- merge(cti_PA_FIN, sites_FIN, by.x = "Site", by.y = "Site")

cti_PA_butterflies <- as.tbl(bind_rows(cbind.data.frame(cti_PA_NL, country = "NL"), cbind.data.frame(cti_PA_FIN, country = "FIN")))
cti_PA_butterflies$Site <- paste0(cti_PA_butterflies$Site, "_", cti_PA_butterflies$country)

# add grid cell (50x50km)
but.pts <- SpatialPoints(cti_AB_butterflies[,c("X", "Y")])
but.r <- raster(ext = extent(but.pts)*1.1, resolution = 50000)
values(but.r) <- c(1:ncell(but.r))
gridCell50 <- extract(but.r, but.pts)

cti_AB_butterflies <- bind_cols(cti_AB_butterflies, gridCell50 = gridCell50)
cti_PA_butterflies <- bind_cols(cti_PA_butterflies, gridCell50 = gridCell50)

cti_butterflies <- as.tbl(bind_rows(cbind.data.frame(type = "Presence", cti_PA_butterflies), cbind.data.frame(type = "Abundance", cti_AB_butterflies)))

# add fragmentation
cti_butterflies <- left_join(cti_butterflies %>% ungroup(), frag_data) 

cti_butterflies <- cti_butterflies %>% group_by(coords = paste(X, Y), Year, type, Scale) %>% summarise_all(first)

cti_butterflies %>% group_by(Site, Year, type, Scale) %>% summarise(n = n()) %>% ungroup() %>% summarise(max(n))
cti_butterflies %>% group_by(coords, Year, type, Scale) %>% summarise(n = n()) %>% ungroup() %>% summarise(max(n))

cti_butterflies <- cti_butterflies %>% ungroup()  %>% dplyr::select(-coords)

# merge and write
write_csv(cti_butterflies, "../Data/cti_butterflies_data_generalistLanduse.csv")



#############
### birds ###
#############

# sites data
sites_SWE <- read_csv(file = "../Data/Birds - Sweden/Sites_SWE_ETRS89_landcover.csv")

# fragmentation data
frag_data_SWE <- read_csv("../Connectivity/Fragmentation/SWE/Frag_indices.csv")

# temperature data
temp_data_SWE <- read_csv("../Data/temperature_SWE.csv")

#dbMEM data
# dbMEM_data_SWE <- read_csv("../Data/Birds - Sweden/dbMEM.SWE.csv")

sites_SWE <- left_join(sites_SWE, frag_data_SWE, by = c("Site" = "LID"))
sites_SWE <- left_join(sites_SWE, lc_class, by = c("Landcover" = "GRID_CODE"))
sites_SWE <- left_join(sites_SWE, temp_data_SWE, by = c("Site"))
# sites_SWE <- left_join(sites_SWE, dbMEM_data_SWE)
sites_SWE <- sites_SWE %>% dplyr:::select(-CLC_CODE,-RGB)

# cti data
# Abundance-weighted #
cti_AB_SWE <- read.csv("../Data/Birds - Sweden/CTI_abundance_Sweden_1996-2016.csv", dec = ".")
cti_AB_SWE <- as.tbl(inner_join(cti_AB_SWE, sites_SWE, by = "Site"))

# P/A #
cti_PA_SWE <- read.csv("../Data/Birds - Sweden/CTI_presence_Sweden_1996-2016.csv", dec = ".")
cti_PA_SWE <- as.tbl(inner_join(cti_PA_SWE, sites_SWE, by = "Site"))

# add grid cell (50x50km)
bir.pts <- SpatialPoints(cti_AB_SWE[,c("X", "Y")])
bir.r <- raster(ext = extent(bir.pts)*1.1, resolution = 50000)
values(bir.r) <- c(1:ncell(bir.r))
gridCell50 <- extract(bir.r, bir.pts)

cti_AB_SWE <- bind_cols(cti_AB_SWE, gridCell50 = gridCell50)
cti_PA_SWE <- bind_cols(cti_PA_SWE, gridCell50 = gridCell50)

# merge and write
cti_birds <- as.tbl(bind_rows(cbind.data.frame(type = "Presence", cti_PA_SWE), cbind.data.frame(type = "Abundance", cti_AB_SWE)))
write_csv(cti_birds, "../Data/cti_birds_data.csv")


###############################################
###############################################
######                                  #######
######      Presence-absence data       #######
######                                  #######
###############################################
###############################################


library(data.table)
library(raster)

#birds
countsSWE1 <- rio::import("../Data/Birds - Sweden/public_totalstandard_Ia.xlsx")
countsSWE2 <- rio::import("../Data/Birds - Sweden/public_totalstandard_Ib.xlsx")
countsSWE3 <- rio::import("../Data/Birds - Sweden/public_totalstandard_IIa.xlsx")
countsSWE4 <- rio::import("../Data/Birds - Sweden/public_totalstandard_IIb.xlsx")
countsSWE <- as.tbl(bind_rows(countsSWE1, countsSWE2, countsSWE3, countsSWE4))

rm(list = c("countsSWE1", "countsSWE2", "countsSWE3", "countsSWE4"))

countsSWE$art <- as.numeric(countsSWE$art)
countsSWE <- countsSWE %>% dplyr:::filter(art %in% 1:645)
countsSWE <- countsSWE[,c(2,4,5,22,23)]

countsSWE$n <- countsSWE$pkind + countsSWE$lind
countsSWE <- countsSWE[,-c(4,5)]

colnames(countsSWE) <- gsub ("karta", "Site", colnames(countsSWE))
colnames(countsSWE) <- gsub ("art", "Species", colnames(countsSWE))
colnames(countsSWE) <- gsub ("yr", "Year", colnames(countsSWE))


countsSWE <- countsSWE %>%
  group_by(Site, Species, Year) %>%
  summarise(n = mean(n)) %>%
  ungroup() %>%
  group_by(Site) %>%
  complete(Year, nesting(Site, Species)) %>%
  mutate(n = ifelse(is.na(n), 0, 1))


#butterflies - finland
countsFIN1 <- fread("../Data/Butterflies - Finland/FINLAND_Records_1999-2015.txt", sep = ";", h=T)
countsFIN2 <- fread("../Data/Butterflies - Finland/FINLAND_Records_2016.txt", sep = ";", h=T)[,-6]
countsFIN <- rbind(countsFIN1, countsFIN2)
colnames(countsFIN) <- gsub ("Species_Faunaeur", "Species", colnames(countsFIN))
colnames(countsFIN) <- gsub ("Individuals", "n", colnames(countsFIN))

countsFIN$Species <- gsub("Polygonia c-album", "Nymphalis c-album", countsFIN$Species)
countsFIN$Species <- gsub("Cyaniris semiargus", "Polyommatus semiargus", countsFIN$Species)
countsFIN$Species <- gsub("Leptidea juvernica", "Leptidea sinapis", countsFIN$Species)


countsFIN <- countsFIN %>% 
  group_by(Site, Species, Year) %>% 
  summarise(n = mean(n)) %>% 
  ungroup() %>%
  group_by(Site) %>%
  complete(Year, nesting(Site, Species)) %>%
  mutate(n = ifelse(is.na(n), 0, 1))

table((countsFIN %>% 
         group_by(Site, Species) %>% summarise(n = length(unique(Year))))$n)


#butterflies - netherlands
# occupancy_NL <- read_csv2("../Data/Butterflies - Netherlands/Occupancy_km_Spec_Yr.csv")

countsNL1 <- readRDS("../Data/Butterflies - Netherlands/AllSpecies_reg_gam_ind_20171206_algroutes.rds") %>% as.tbl %>%
  dplyr::select(1,2,3,4) %>% rename(n = regional_gam)
countsNL2 <- rio::import("../Data/Butterflies - Netherlands/MissingSpecies.xlsx") %>% as.tbl %>%
  dplyr::select(1,3,4,5) %>% rename(n = Ntot, SITE = Site) %>% mutate(n = ifelse(n == -1, 0 , n))
countsNL <- bind_rows(countsNL1, countsNL2)

countsNL <- countsNL %>% 
  group_by(SITE) %>%
  complete(YEAR, nesting(SITE, SPECIES)) %>%
  mutate(n = ifelse(is.na(n) | n == 0, 0, 1)) %>%
  rename(Year = YEAR, Site = SITE, Species = SPECIES) %>% dplyr::select(Year,Site,Species,n)

table((countsNL %>% 
         group_by(Site, Species) %>% summarise(n = length(unique(Year))))$n)

#################################
#### merge and write on disc ####
#################################

#birds
# birds_sti <- read.csv("../Birds Atlas/BOTW/STI.csv")
# list.spSWE <- read.csv("../Data/Birds - Sweden/Species_list.csv")
# birds_sti <- merge(birds_sti, list.spSWE, by.x = "Species", by.y = "AltName")[,1:3]
# 
# birds.data <- as.tbl(fread("../Data/cti_birds_data.csv"))
# birds.data <- birds.data %>% dplyr:::filter(Scale %in% c(3000,30000), type == "Presence")
# birds.data <- birds.data %>% group_by(Site, Scale) %>% dplyr:::filter(row_number(type) == 1)
# birds.data <- birds.data[,c(3,5,6,8,10,11,17,18)]
# 
# pres_abs.data <- left_join(countsSWE, birds.data, by = c("Site" = "Site"))
# pres_abs.data <- left_join(pres_abs.data, birds_sti, by = c("Species" = "Species")) %>% mutate(Species = Species.y) %>% dplyr:::select(-Species.y)
# 
# write_csv(pres_abs.data, "../Data/pres_abs_birds_data.csv")

frag_data <- read_csv("../Connectivity - Fragmentation/Fragmentation/Frag_indices_Allhab.csv")

traits.butterflies <- as.tbl(read_csv("../Data/Traits/SpeciesTraits_WDV2014.csv"))
traits.butterflies[grep("Colias croceus", traits.butterflies$Scientific_name), "Scientific_name"] <- "Colias crocea"
traits.butterflies[grep("walbum", traits.butterflies$Scientific_name), "Scientific_name"] <- "Satyrium w-album"
traits.butterflies[grep("Neozephyrus quercus", traits.butterflies$Scientific_name), "Scientific_name"] <- "Favonius quercus"
traits.butterflies[grep("lycaon", traits.butterflies$Scientific_name), "Scientific_name"] <- "Hyponephele lycaon"
traits.butterflies[grep("Polygonia", traits.butterflies$Scientific_name), "Scientific_name"] <- "Nymphalis c-album"
traits.butterflies[grep("Inachis io", traits.butterflies$Scientific_name), "Scientific_name"] <- "Aglais io"
traits.butterflies[grep("tithonus", traits.butterflies$Scientific_name), "Scientific_name"] <- "Pyronia tithonus"
traits.butterflies[grep("alcon", traits.butterflies$Scientific_name), "Scientific_name"] <- "Phengaris alcon"

butterfly_habitat <- read_csv("../Data/Butterfly_biotopes.csv")
butterfly_habitat[grep("Neozephyrus quercus", butterfly_habitat$Species), "Species"] <- "Favonius quercus"

#butterflies - finland

sites_FIN <- read_csv(file = "../Data/Butterflies - Finland/Sites_FIN_ETRS89_landcover.csv")

# add grid cell (50x50km)
but.pts <- SpatialPoints(sites_FIN[,c("X", "Y")])
but.r <- raster(ext = extent(but.pts)*1.1, resolution = 50000)
values(but.r) <- c(1:ncell(but.r))
gridCell50 <- raster::extract(but.r, but.pts)
sites_FIN <- bind_cols(sites_FIN, gridCell50 = gridCell50)

countsFIN <- left_join(countsFIN, sites_FIN, by = c("Site")) %>% ungroup() %>%
  mutate(Site = paste0(Site, "_FIN"))


# FIN.data <- as.tbl(fread("../Data/cti_butterflies_data.csv")) %>% dplyr:::filter(country == "FIN", type == "Presence")
# FIN.data <- FIN.data %>% group_by(Site, Scale) %>% dplyr:::filter(row_number(type) == 1)
# FIN.data <- FIN.data[,-c(3,4)]

pres_abs.data <- left_join(countsFIN,
                           traits.butterflies %>% dplyr:::select(-Species),
                           by = c("Species" = "Scientific_name"))
pres_abs.data <- left_join(pres_abs.data, butterfly_habitat, by = "Species")

# pres_abs.data <- left_join(pres_abs.data, frag_data %>% dplyr::filter(Habitat == "Generalist"), 
#   by = c("Site" = "Site"))

pres_abs.data <- bind_rows(
  left_join(
    pres_abs.data %>% dplyr::filter(Habitat == "Forest"), frag_data,
    by = c("Site" = "Site", "Habitat" = "Habitat")),
  left_join(
    pres_abs.data %>% dplyr::filter(Habitat == "Grassland" | Habitat == "Wetland") %>%
      mutate(Habitat = "Open"),
    frag_data,
    by = c("Site" = "Site", "Habitat" = "Habitat")),
  left_join(
    pres_abs.data %>% dplyr::filter(is.na(Habitat)) %>%
      mutate(Habitat = "Generalist"),
    frag_data,
    by = c("Site" = "Site", "Habitat" = "Habitat")))

pres_abs.data %>% dplyr::filter(is.na(PC1)) %>% ungroup () %>% group_by(Species) %>% summarise(n = length(unique(Site)))

write_csv(pres_abs.data, "../Data/pres_abs_FIN_data.csv")

#butterflies - netherlands
sites_NL <- as.tbl(left_join(rio:::import(file = "../Data/Butterflies - Netherlands/Sites_NL.xlsx"),
                             read_csv(file = "../Data/Butterflies - Netherlands/Sites_NL_ETRS89_landcover.csv")))
# add grid cell (50x50km)
but.pts <- SpatialPoints(sites_NL[,c("X", "Y")])
but.r <- raster(ext = extent(but.pts)*1.1, resolution = 50000)
values(but.r) <- c(1:ncell(but.r))
gridCell50 <- raster::extract(but.r, but.pts)

sites_NL <- bind_cols(sites_NL, gridCell50 = gridCell50)

countsNL <- left_join(countsNL, sites_NL %>% dplyr::select(-X_km, -Y_km)) %>% ungroup() %>%
  mutate(Site = paste0(Site, "_NL")) %>% filter(!is.na(X))

pres_abs.data <- left_join(countsNL,
                           traits.butterflies %>% mutate(Species = casefold(Species))) %>%
  dplyr::mutate(Species = Scientific_name) %>% dplyr::select(-Scientific_name) %>% filter(!is.na(Species))

pres_abs.data <- left_join(pres_abs.data, butterfly_habitat, by = "Species")

# pres_abs.data <- left_join(pres_abs.data, frag_data %>% dplyr::filter(Habitat == "Generalist"), 
#                            by = c("Site" = "Site"))

pres_abs.data <- bind_rows(
  left_join(
    pres_abs.data %>% dplyr::filter(Habitat == "Forest"), frag_data,
    by = c("Site" = "Site", "Habitat" = "Habitat")),
  left_join(
    pres_abs.data %>% dplyr::filter(Habitat == "Grassland" | Habitat == "Wetland") %>%
      mutate(Habitat = "Open"),
    frag_data,
    by = c("Site" = "Site", "Habitat" = "Habitat")),
  left_join(
    pres_abs.data %>% dplyr::filter(is.na(Habitat)) %>%
      mutate(Habitat = "Generalist"),
    frag_data,
    by = c("Site" = "Site", "Habitat" = "Habitat")))

write_csv(pres_abs.data, "../Data/pres_abs_NL_data.csv")



#### add relative STI

butterflies.cti.presence <- as.tbl(fread("../Data/cti_butterflies_data.csv")) %>% dplyr::filter(type == "Presence")
butterflies <- bind_rows(read_csv("../Data/pres_abs_FIN_data.csv"), read_csv("../Data/pres_abs_NL_data.csv"))

butterflies <- left_join(butterflies,
                         butterflies.cti.presence %>% group_by(Site) %>% summarise(cti = mean(cti)))

butterflies <- butterflies %>% mutate(STI_rel = cti - STI) %>%
  dplyr::filter(!is.na(STI)) %>%
  mutate(gridCell50 = ifelse(country == "NL", paste0(gridCell50, "_NL"), paste0(gridCell50, "_FIN")))

write_csv(butterflies, "C:/Local Folder (c)/butterflies_occ.csv")
