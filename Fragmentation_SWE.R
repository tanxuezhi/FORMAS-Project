library(raster)
library(tidyverse)
library(rgeos)

birds.data <- read_csv("../Data/Birds - Sweden/Sites_SWE_ETRS89_landcover.csv")

# load swedish birds monitoring sites
sites.SWE <- rio::import("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Data/Birds - Sweden/Koordinater standardrutternas samtliga punkter.xlsx")
sites.SWE <- sites.SWE %>% group_by(karta) %>% summarise(x = mean(rt90_o), y = mean(rt90_n)) 

# load CLC
CLC_SWE <- raster("../Landcover/Corine_land-cover_2012_raster/CLC_SWE.tif")

# create mask for monitoring sites
r <- sites.SWE %>% complete(x, y) %>% dplyr::select(-karta)
r <- SpatialPoints(r)
crs(r) <- "+init=epsg:3021"
r <- r@coords
r <- rasterFromXYZ(cbind.data.frame(r,y = 1:nrow(r)))

mask <- setdiff(Which(r > 0, cells = TRUE), raster::extract(r, sites.SWE[,2:3], cellnumbers = T)[,"cells"])
r[mask] <- NA

grid <- rasterToPolygons(r)
crs(grid) <- "+init=epsg:3021"
grid <- spTransform(grid, crs(CLC_SWE))


# reclass
LC_reclass_open <- as.tbl(rio::import(file = "../Landcover/Corine_land-cover_2012_raster/clc_legend.xls", which = 5L))
LC_reclass_forest<- as.tbl(rio::import(file = "../Landcover/Corine_land-cover_2012_raster/clc_legend.xls", which = 6L))

LC_reclass <- cbind.data.frame(LC_reclass_open[,1:3], LC_reclass_open[,4] + LC_reclass_forest[,4])

SNH <- reclassify(CLC_SWE, cbind.data.frame(is = LC_reclass[,1], becomes = LC_reclass[,4]))
writeRaster(SNH, "../Landcover/SNH/SNH_SWE.tif")


#### fragmentation
frag <- c()
for(i in 1:length(grid)){
  
  print(i / length(grid) * 100)
  
  r <- crop(SNH, grid[i,])

  ai = length(Which(r == 1, cells = TRUE))
  
  if(length(r[!is.na(r)]) - ai > 1 & ai > 0){
    
    a <- adjacent(r, 1:ncell(r), 4, pairs=TRUE)
    tb <- table(r[a[,1]], r[a[,2]])
    
    gii = tb[2,2]
    gik = sum(tb[2,])
    
    P = ai/length(r[!is.na(r)])
    
    G = gii / gik
    
    CLUMPY = ifelse(G < P & P < .5,(G - P)/P, (G - P)/(1 - P))
  }
  
  frag <- rbind.data.frame(frag, cbind.data.frame(gCentroid(grid[i,])@coords, PLAND = P, CLUMPY))

}

# merge to monitoring sites
grid@data <- frag
grid <- spTransform(grid, crs("+init=epsg:3021"))

frag_data <- cbind.data.frame(sites.SWE, raster::intersect(SpatialPoints(sites.SWE[,-1]), grid)@data[,3:4])

frag_raster <- rasterFromXYZ(frag_data[,-1])


# write
write_csv(frag_data, "../Connectivity - Fragmentation/Fragmentation/frag_SWE.csv")
writeRaster(frag_raster[[1]], "../Connectivity - Fragmentation/Fragmentation/PLAND_SWE.tif")
writeRaster(frag_raster[[2]], "../Connectivity - Fragmentation/Fragmentation/CLUMPY_SWE.tif")

