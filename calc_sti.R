library(dismo)
library(raster)
library(sp)
library(rgdal)
library(data.table)
library(dplyr)
library(rgeos)

#### temperature data ####
# load temperature data
clim.Folder <- "//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Climate data/"
breeding.Months <- c(3:8)
europe.Extent <- extent(c(-32,70,30,82))

temperature <- getData("worldclim", var='tmean', res = 5, path = clim.Folder)
temperature <- crop(x=temperature, y=europe.Extent)
temperature.breeding <- temperature[[breeding.Months]]
temperature.breeding <- mean(temperature.breeding)/10

writeRaster(temperature.breeding, paste0(clim.Folder, "/MeanTemp_March_august.tif"), overwrite = T)


#### Load data ####
temperature.breeding <- raster("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Climate data/MeanTemp_March_august.tif")
temperature.breeding <- aggregate(temperature.breeding, 4)

# load functions
source("functions.R")

# load species lists (monitoring and Birdlife shapefiles)
species_BOTW <- read.csv("../Data/Birds Atlas/BOTW/Species.csv")[,1:6]
species_EBCC <- read.csv("../Data/Birds Atlas/EBCC1/species_EBCC1.csv")

species <- read.csv(file = "../sp.csv")

species[!species$AltName %in% species_BOTW$ScientificName, "AltName"]
species[!species[,1] %in% species_EBCC$Species,]

species %>% filter(!AltName %in% species_EBCC$Species & !latin %in% species_EBCC$Species) %>% select(3)

syn <- synonyms(species[!species[,1] %in% species_EBCC$Species,][1:2], db = "iucn", rows = 1)

#### Using BOTW ####

BOTW <- readOGR(dsn = "../Data/Birds Atlas/BOTW/Selected_Species_breedORresid_dis.shp")

temp_extract <- extract(temperature.breeding, BOTW)

sti <- bind_cols(Species = BOTW$SCINAME, sti = (unlist(lapply(temp_extract, function(x)mean(x, na.rm = T))))) %>%
  group_by(Species) %>% summarise(sti = mean(sti, na.rm = T))

write.csv(sti, "../Data/Birds Atlas/BOTW/STI.csv", row.names = F)

#### Using EBCC1 ####

EBCC <- fread("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Data/Birds Atlas/EBCC1/0005875-170627171947987.csv")

EBCC$species <- gsub("Tetrastes bonasia", "Bonasa bonasia", EBCC$species)
EBCC$species <- gsub("Linaria cannabina", "Carduelis cannabina", EBCC$species)
EBCC$species <- gsub("Acanthis flammea", "Carduelis flammea", EBCC$species)
EBCC$species <- gsub("Linaria flavirostris", "Carduelis flavirostris", EBCC$species)
EBCC$species <- gsub("Spinus spinus", "Carduelis spinus", EBCC$species)
EBCC$species <- gsub("Corvus corone", "Corvus corone cornix", EBCC$species)
EBCC$species <- gsub("Coloeus monedula", "Corvus monedula", EBCC$species)
EBCC$species <- gsub("Dryobates minor", "Dendrocopos leucotos", EBCC$species)
EBCC$species <- gsub("Capella media", "Gallinago media", EBCC$species)
EBCC$species <- gsub("Lagopus muta", "Lagopus mutus", EBCC$species)
EBCC$species <- gsub("Mergellus albellus", "Mergus albellus", EBCC$species)
EBCC$species <- gsub("Poecile cinctus", "Parus cinctus", EBCC$species)
EBCC$species <- gsub("Poecile montanus", "Parus montanus", EBCC$species)
EBCC$species <- gsub("Poecile palustris", "Parus palustris", EBCC$species)
EBCC$species <- gsub("Regulus ignicapilla", "Regulus ignicapillus", EBCC$species)
EBCC$species <- gsub("Sternula albifrons", "Sterna albifrons", EBCC$species)
EBCC$species <- gsub("Hydroprogne caspia", "Sterna caspia", EBCC$species)
EBCC$species <- gsub("Thalasseus sandvicensis", "Sterna sandvicensis", EBCC$species)
EBCC$species <- gsub("Morus bassanus", "Sula bassana", EBCC$species)
EBCC$species <- gsub("Lyrurus tetrix", "Tetrao tetrix", EBCC$species)
EBCC$species <- gsub("Acanthis hornemanni", "Carduelis hornemanni", EBCC$species)
EBCC$species <- gsub("Picoides leucotos", "Dendrocopos leucotos", EBCC$species)
EBCC$species <- gsub("Larus graellsi", "Larus fuscus", EBCC$species)
EBCC$species <- gsub("Muscicapa striata", "Muscicapa striata", EBCC$species)
EBCC[EBCC$scientificname == "Stictocarbo Bonaparte, 1855","species"] <- "Phalacrocorax aristotelis"
EBCC$species <- gsub("Phylloscopus sibillatrix", "Phylloscopus sibilatrix", EBCC$species)
EBCC$species <- gsub("Setophaga striata", "Muscicapa striata", EBCC$species)
EBCC[EBCC$genus == "Mareca","species"] <- "Anas strepera"
EBCC[EBCC$genus == "Erithacus","species"] <- "Erithacus rubecula"
EBCC <- rbind(EBCC, EBCC[EBCC$genus == "Mareca",] %>% mutate(species = "Anas penelope"))
EBCC$species <- gsub("Hydrocoloeus minutus", "Larus minutus", EBCC$species)
EBCC$species <- gsub("Chroicocephalus ridibundus", "Larus ridibundus", EBCC$species)


sort(species[!species$species %in% EBCC$species, ])

sti.list <- EBCC %>% group_by(species) %>% 
  do(sti = as.data.frame(mean(raster::extract(temperature.breeding, 
                                              .[,c("decimallongitude","decimallatitude")]), na.rm =T))) %>%
  unnest
colnames(sti.list)[2] <- "sti"
sti.list <- sti.list %>% filter(species != "") 
sti.list <- bind_rows(sti.list, data.frame(species = "Loxia species", 
                                 sti = sti.list %>% filter(grepl("Loxia", species)) %>% summarise(sti = mean(sti))))

sti.list.sd <- EBCC %>% group_by(species) %>% 
  do(sti = as.data.frame(sd(raster::extract(temperature.breeding, 
                                              .[,c("decimallongitude","decimallatitude")]), na.rm =T))) %>%
  unnest
colnames(sti.list.sd)[2] <- "sti.sd"
sti.list.sd <- sti.list.sd %>% filter(species != "") 
sti.list.sd <- bind_rows(sti.list.sd, data.frame(species = "Loxia species", 
                                           sti = sti.list.sd %>% filter(grepl("Loxia", species)) %>% 
                                             summarise(sti.sd = mean(sti.sd, na.rm =T))))
sti.list.sd <- sti.list.sd %>% mutate(sti.sd = ifelse(is.na(sti.sd), 0, sti.sd))

sti.list <- sti.list %>% left_join(sti.list.sd)

write_csv(sti.list, "//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Data/Birds Atlas/EBCC1/EBCC1_sti.csv")

#### using IUCN data ####
IUCN_distri <- readOGR("../Birds Atlas/BOTW/Selected_Species_breedORresid_dis.shp")
# IUCN_distri <- crop(IUCN_distri, extent(temperature.breeding))
save(IUCN_distri, file = "IUCN_distri.RData")
load("IUCN_distri.RData")


#loop over species and calculate Species Temperature Index
sti.list <- sti(distri = IUCN_distri, temperature = temperature.breeding, output = "../Birds Atlas/BOTW/STI.csv")


##### using GBIF #####
library(rgbif)


occ_search(scientificName = as.character(species[1,]), hasCoordinate = T, hasGeospatialIssue = F,
           decimalLongitude = paste(extent(temperature.breeding)@xmin, ",", extent(temperature.breeding)@xmax),
           decimalLatitude = paste(extent(temperature.breeding)@ymin, ",", extent(temperature.breeding)@ymax),
           geometry = wkt)
