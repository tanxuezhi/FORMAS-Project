library(dismo)
library(raster)
library(sp)
library(rgdal)

# load functions
source("functions.R")

# load species lists (monitoring and Birdlife shapefiles)
species <- rio::import(file = "//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Data/Birds - Sweden/public_eurolist.xlsx", which = 1L)
IUCN_distri <- list.files("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/IUCN data/Birdlife distributions", 
                          pattern = ".shp", 
                          recursive = T,
                          full.names = T)

# load temperature data
clim.Folder <- "//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Climate data/"
breeding.Months <- c(3:8)
europe.Extent <- extent(c(-15,50,35,75))

temperature <- getData("worldclim", var='tmean', res = 5, path = clim.Folder)
temperature <- crop(x=temperature, y=europe.Extent)
temperature.breeding <- temperature[[breeding.Months]]
temperature.breeding <- mean(temperature.breeding)/10

# species names
species.names <- species[,"latin"]
species.names <- gsub(" ", "_", species.names)

# select seasons
# 1 = resident
# 2 = breeding
season <- c(1,2)

#loop over species and calculate Species Temperature Index
sti.list <- sti(species = species.names, distri = IUCN_distri, temperature = temperature.breeding)
