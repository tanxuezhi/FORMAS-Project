library(dismo)
library(raster)
library(sp)
library(rgdal)
library(data.table)
library(dplyr)

# load functions
source("functions.R")

# load species lists (monitoring and Birdlife shapefiles)

species_distri <- read.csv("../Birds Atlas/BOTW/Species.csv")[,1:6]

species <- read.csv(file = "../Data/Birds - Sweden/Species_list.csv")[1:385,]

nrow(species[species$AltName %in% species_distri$ScientificName, ])
nrow(species)

species[!species$AltName %in% species_distri$ScientificName, "AltName"]

species_distri[grep("Acanthis", species_distri$ScientificName), "ScientificName"]
species_distri[grep("Gadwall", species_distri$CommonName), "ScientificName"]

#### using IUCN data ####
IUCN_distri <- readOGR("../Birds Atlas/BOTW/Selected_Species_breedORresid_dis.shp")
# IUCN_distri <- crop(IUCN_distri, extent(temperature.breeding))
save(IUCN_distri, file = "IUCN_distri.RData")
load("IUCN_distri.RData")

# load temperature data
clim.Folder <- "//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Climate data/"
breeding.Months <- c(3:8)
# europe.Extent <- extent(c(-11,34,32,72))
europe.Extent <- extent(c(-32,70,30,82))

temperature <- getData("worldclim", var='tmean', res = 5, path = clim.Folder)
temperature <- crop(x=temperature, y=europe.Extent)
temperature.breeding <- temperature[[breeding.Months]]
temperature.breeding <- mean(temperature.breeding)/10

writeRaster(temperature.breeding, paste0(clim.Folder, "/MeanTemp_March_august.tif"), overwrite = T)

#loop over species and calculate Species Temperature Index
sti.list <- sti(distri = IUCN_distri, temperature = temperature.breeding, output = "../Birds Atlas/BOTW/STI.csv")


#### Using EBCC1 ####
countsSWE1 <- rio::import("../Data/Birds - Sweden/public_totalstandard_Ia.xlsx")
countsSWE2 <- rio::import("../Data/Birds - Sweden/public_totalstandard_Ib.xlsx")
countsSWE3 <- rio::import("../Data/Birds - Sweden/public_totalstandard_IIa.xlsx")
countsSWE4 <- rio::import("../Data/Birds - Sweden/public_totalstandard_IIb.xlsx")

countsSWE <- as.tbl(bind_rows(countsSWE1, countsSWE2, countsSWE3, countsSWE4))
countsSWE$art <- as.numeric(countsSWE$art)
countsSWE <- countsSWE %>% filter(art %in% 1:645)

species <- read.csv(file = "../Data/Birds - Sweden/Species_list.csv")[1:385,]
species <- species %>% filter(art %in% countsSWE$art)


EBCC <- as.tbl(fread("../EBCC1/0005875-170627171947987.csv")[,c("species", "decimallongitude", "decimallatitude")])
sel.sp <- EBCC %>% filter(species == "")

min(EBCC$decimallongitude, na.rm =T)

sti.list <- EBCC %>% group_by(species) %>% do(as.data.frame(mean(extract(temperature.breeding, .[,2:3]), na.rm =T)))

EBCC %>% filter(., grepl("penelope", species))

species[!species$AltName %in% sti.list$species & !species$latin %in% sti.list$species, "AltName"]
sti.list[grep("Anas", sti.list$species), "species"]
