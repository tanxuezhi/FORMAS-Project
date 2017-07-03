# load libraries
library(raster)
library(SDMTools)
library(rgdal)
library(reshape2)
library(rgeos)

###load sites
sites <- read.csv(file = "../Data/Butterflies - Netherlands/Sites_ETRS89_landcover.csv")
lc_class <- rio::import(file = "../Landcover/clc_legend.xls", which = 1L)
sites <- merge(sites, lc_class, by.x = "Landcover", by.y = "CLC_CODE")

###load landcover
#raster
CLC_NL <- raster("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Landcover/Corine_land-cover_2012_raster/g100_clc12_V18_5.tif")

CLC_NL_SNH <- raster("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Landcover/Corine_land-cover_2012_raster/CLC_2012_100m_reclass_snh_NL.tif")

#load roads
roads <- readOGR("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/EGM_9-0SHP_20161017/DATA/Countries/NL/RoadL.shp")

#shapefile
CLC_NL <- readOGR("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Landcover/CLC_2012_NL.shp")

#landuse stats at monitoring sites
table(sites$LABEL2)

######## prepare inputs for UNICOR ########
# set directory
anal_dir <- "C:/Local Folder (c)/Connectivity softwares/UNICOR-master/Real_test_NL"

reclassTable <- rio::import(file = "../Landcover/Corine_land-cover_2012/Legend/clc_legend.xls", which = 3L)
resistanceForest <- reclassify(CLC_NL, reclassTable[,-1])
resistanceForest <- crop(resistanceForest, extent(CLC_NL)/5)
writeRaster(resistanceForest, paste0(anal_dir,"/NL_resistance_forest.asc"), overwrite=T)

resistanceForestTxt  <- readLines(paste0(anal_dir,"/NL_resistance_forest.asc"))
resistanceForestTxt  <- gsub(pattern = "NROWS", replace = "nrows", x = resistanceForestTxt)
resistanceForestTxt  <- gsub(pattern = "NCOLS", replace = "ncols", x = resistanceForestTxt)
resistanceForestTxt  <- gsub(pattern = "XLLCORNER", replace = "xllcorner", x = resistanceForestTxt)
resistanceForestTxt  <- gsub(pattern = "YLLCORNER", replace = "yllcorner", x = resistanceForestTxt)
resistanceForestTxt  <- gsub(pattern = "CELLSIZE", replace = "cellsize", x = resistanceForestTxt)
writeLines(resistanceForestTxt, con=paste0(anal_dir,"/NL_resistance_forest.asc"))


# sample points in selected LC
Forest <- resistanceForest
Forest[!Forest %in% 1] <- NA

Forest <- rasterToPolygons(Forest, dissolve = T)
Forest <- disaggregate(Forest)

pts <- gCentroid(Forest, byid = T)
write.csv(pts@coords, paste0(anal_dir,"/samplePtsForest.xy"), row.names = F)

# euclidean distance
dist_btw_patches <- as.matrix(dist(pts@coords))
dist_btw_patches[upper.tri(dist_btw_patches)] <- NA
dist_btw_patches <- melt(as.matrix(dist_btw_patches), varnames = c("pt1", "pt2"), value.name = "Distance")
dist_btw_patches <- na.omit(dist_btw_patches)
dist_btw_patches[dist_btw_patches[,3]==0,3] <- NA
dist_btw_patches <- na.omit(dist_btw_patches)

write.table(dist_btw_patches, paste0(anal_dir, "/eucl_distance.txt"), row.names = F, col.names = F)


## run UNICOR
#python path
py <- "C:/Anaconda/python.exe"
#UNICOR path
unicor <- "C:/Local Folder (c)/Connectivity softwares/UNICOR-master/unicor/unicor.py"

backWd <- getwd()
setwd(anal_dir)
system(paste0(py, " \"", unicor, "\" ", " \"", anal_dir, "/Real_rest_NL.rip", "\""))
setwd(backWd)

## result
paths <- raster(paste0(anal_dir,"/NL_resistance_forest.asc_samplePtsForest.addedpaths.txt"))

dist_btw_patches <- read.csv(paste0(anal_dir,"/NL_resistance_forest.asc_samplePtsForest.cdmatrix.csv"), h = F)
colnames(dist_btw_patches) <- gsub("V", "", colnames(dist_btw_patches))
dist_btw_patches[upper.tri(dist_btw_patches)] <- NA
dist_btw_patches <- melt(as.matrix(dist_btw_patches), varnames = c("pt1", "pt2"), value.name = "Distance")
dist_btw_patches <- na.omit(dist_btw_patches)
dist_btw_patches[dist_btw_patches[,3]==0,3] <- NA
dist_btw_patches <- na.omit(dist_btw_patches)

write.table(dist_btw_patches, paste0(anal_dir, "/distance.txt"), row.names = F, col.names = F)
write.table(cbind(1:length(Forest), as.data.frame(gArea(Forest, byid=T))), paste0(anal_dir, "/nodes.txt"), row.names = F, col.names = F)

## run conefor
# find minumum threshold
max(aggregate(Distance ~ pt1, dist_btw_patches, min))

#conefor path
conefor <- "C:/Local Folder (c)/Connectivity softwares/Conefor/Conefor_command_line/Conefor_2014_Windows_32_and_64_bit/coneforWin64.exe"
  
backWd <- getwd()
setwd(anal_dir)
system(paste0("\"", conefor,"\"", " -nodeFile ", " \"", anal_dir, "/nodes.txt", "\"", 
              " -conFile ", " \"", anal_dir, "/eucl_distance.txt", "\" -t dist -confAdj 10000 -BC"))
setwd(backWd)

## result2
nodeImpEucl <- read.table(paste0(anal_dir, "/node_importances_Eucl_10000.txt"), h=T)
nodeImpLCP <- read.table(paste0(anal_dir, "/node_importances_LCP_50000.txt"), h=T)

par(mfrow=c(1,2))
plot(nodeImpEucl$BC ~ nodeImpLCP$BC, main = "BC", ylab = "Euclidean distance", xlab = "Least-cost path distance")
plot(nodeImpEucl$dIIC ~ nodeImpLCP$dIIC, main = "dIIC", ylab = "Euclidean distance", xlab = "Least-cost path distance")
par(mfrow=c(1,1))


#plot
par(mfrow=c(1,2))
plot(paths, legend = F, axes = F, box = F)
plot(Forest, add=T, border=NA,
     col = colorRampPalette(c('blue','red'))(length(pts))[as.numeric(cut((nodeImpEucl$dIIC^(1/2)),breaks = length(pts)))])
plot(paths, legend = F, axes = F, box = F)
plot(Forest, add=T, border=NA,
     col = colorRampPalette(c('blue','red'))(length(pts))[as.numeric(cut((nodeImpEucl$BC^(1/2)),breaks = length(pts)))])
par(mfrow=c(1,1))

