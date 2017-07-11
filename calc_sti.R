library(dismo)
library(raster)
library(sp)
library(rgdal)

# load functions
source("functions.R")

# load species lists (monitoring and Birdlife shapefiles)

species_distri <- read.csv("../BOTW/Species.csv")[,1:6]
species <- read.csv(file = "../Data/Birds - Sweden/Species_list.csv")[1:385,]

nrow(species[species$AltName %in% species_distri$ScientificName, ])
nrow(species)

species[!species$AltName %in% species_distri$ScientificName, "AltName"]

species_distri[grep("Acanthis", species_distri$ScientificName), "ScientificName"]
species_distri[grep("Gadwall", species_distri$CommonName), "ScientificName"]

IUCN_distri <- readOGR("../BOTW/Selected_Species_breedORresid_dis.shp")
# IUCN_distri <- crop(IUCN_distri, extent(temperature.breeding))
save(IUCN_distri, file = "IUCN_distri.RData")
load("IUCN_distri.RData")

# load temperature data
clim.Folder <- "//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Climate data/"
breeding.Months <- c(3:8)
europe.Extent <- extent(c(-11,34,32,72))

temperature <- getData("worldclim", var='tmean', res = 5, path = clim.Folder)
temperature <- crop(x=temperature, y=europe.Extent)
temperature.breeding <- temperature[[breeding.Months]]
temperature.breeding <- mean(temperature.breeding)/10

writeRaster(temperature.breeding, paste0(clim.Folder, "/MeanTemp_March_august.tif"))

#loop over species and calculate Species Temperature Index
sti.list <- sti(distri = IUCN_distri, temperature = temperature.breeding, output = "../BOTW/STI.csv")
