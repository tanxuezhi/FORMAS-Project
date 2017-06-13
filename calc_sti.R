library(dismo)
library(rgdal)

wd <- "//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Data/Birds - Sweden"
setwd(wd)

# load species lists (monitoring and Birdlife shapefiles)
species <- rio::import(file = "//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Data/Birds - Sweden/public_eurolist.xlsx", which = 1L)
IUCN_distri <- list.files("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/IUCN data/Birdlife distributions", 
                          pattern = ".shp", 
                          recursive = T,
                          full.names = T)

# load temperature data
temperature <- getData("worldclim", var='tmean', res = 5)
temperature <- crop(x=temperature, y=extent(c(-15,50,35,75)))
temperature.breeding <- temperature[[c(3:8)]]
temperature.breeding <- mean(temperature.breeding)/10

# species names
species.names <- species[,"latin"]
species.names <- gsub(" ", "_", species.names)

# select seasons
# 1 = resident
# 2 = breeding
season <- c(1,2)

#loop over species and calculate Species Temperature Index
sti_list <- c()
for(i in species.names){
  cat(paste("Species:", i, "\n"))
  nb.sp <- grep(i, IUCN_distri)
  
  if(length(nb.sp)>0){
  cat(paste("Found in birdlife data", "\n"))
  distri_shp <- readOGR(IUCN_distri[nb.sp][1], verbose = F) # load shapefile
  distri_shp <- distri_shp[distri_shp@data$SEASONAL %in% season,] # select resident and breeding ranges
  
  sti <- mean(unlist(extract(temperature.breeding, SpatialPolygons(distri_shp@polygons))), na.rm=T) # calculate STI
  cat(paste("STI =", round(sti, 2), "?C", "\n\n"))
  sti_list <- rbind.data.frame(sti_list, cbind.data.frame(Species = i, STI = sti))
  }else{
    cat(paste("Not found in birdlife data", "\n\n"))
    sti_list <- rbind.data.frame(sti_list, cbind.data.frame(Species = i, STI = NA))
  }
}


