library(corrplot)
library(raster)
library(sp)
library(stringr)
library(rgeos)

###load landcover
# reclassified by semi-natural habitat
CLC_SNH <- raster("../Landcover/SNH/SNH_merged.tif")

###load sites and write cropped SNH rasters
sites_FIN <- read.csv(file = "../Data/Butterflies - Finland/Sites_FIN_ETRS89_landcover.csv")
sites_FIN <- SpatialPointsDataFrame(sites_FIN[,2:3], sites_FIN, proj4string = CRS("+init=epsg:3035"))
CLC_SNH_FIN <- crop(CLC_SNH, extent(sites_FIN) + 100000)
writeRaster(CLC_SNH_FIN, "../Landcover/Finland/SNH_Fin.tif")

sites_NL <- read.csv(file = "../Data/Butterflies - Netherlands/Sites_NL_ETRS89_landcover.csv")
sites_NL <- SpatialPointsDataFrame(sites_NL[,2:3], sites_NL, proj4string = CRS("+init=epsg:3035"))
CLC_SNH_NL <- crop(CLC_SNH, extent(sites_NL) + 100000)
writeRaster(CLC_SNH_NL, "../Landcover/Netherlands/SNH_NL.tif")


##### prepare Fragstat points #####
pts_Fragstat_FIN <- cbind.data.frame(ID=sites_FIN@data$Site, row = rowFromY(CLC_SNH_FIN, sites_FIN@coords[,"Y"]), col = colFromX(CLC_SNH_FIN, sites_FIN@coords[,"X"]))
pts_Fragstat_FIN <- data.frame(paste0("[",pts_Fragstat_FIN[,1],":",pts_Fragstat_FIN[,2],":",pts_Fragstat_FIN[,3],"]"))
colnames(pts_Fragstat_FIN) <- "FPT_TABLE"
write.table(pts_Fragstat_FIN, "../Connectivity/Fragmentation/FIN/FIN_sites.fpt", row.names = F, quote = F)

pts_Fragstat_NL <- cbind.data.frame(ID=sites_NL@data$Site, row = rowFromY(CLC_SNH_NL, sites_NL@coords[,"Y"]), col = colFromX(CLC_SNH_NL, sites_NL@coords[,"X"]))
pts_Fragstat_NL <- data.frame(paste0("[",pts_Fragstat_NL[,1],":",pts_Fragstat_NL[,2],":",pts_Fragstat_NL[,3],"]"))
colnames(pts_Fragstat_NL) <- "FPT_TABLE"
write.table(pts_Fragstat_NL, "../Connectivity/Fragmentation/NL/NL_sites.fpt", row.names = F, quote = F)

##### load results #####
folder <- "../Connectivity/Fragmentation/NL/"
files_frag <- list.files(folder, pattern = ".class", full.names = T)

frag_all_scales <- c()
for(i in 1:length(files_frag)){
  frag <- readLines(files_frag[i])
  frag <- gsub(pattern = "N/A", replace = " ", x = frag)
  writeLines(frag, con=files_frag[i])
  
  frag <- read.table(files_frag[i], h=T, sep = ",", dec=".")
  frag[,1] <- as.numeric(gsub("point_", "", frag[,1]))
  frag <- frag[,c("LID","NLSI","CA", "PLAND")]
  
  frag_all_scales <- rbind.data.frame(frag_all_scales, cbind.data.frame(Scale = str_sub(files_frag[i], -11, -7), frag))
}
frag_all_scales$Scale <- gsub("_","",frag_all_scales$Scale)
write.csv(frag_all_scales, paste0(folder, "Frag_indices.csv"), row.names = F)

# corrplot(abs(cor(frag[,-c(1,2,36)])), method=c("color"),  
#          type="upper")
# plot(frag[,c("NLSI","CA","PROX_MN","CONNECT", "LSI")])
# plot(frag[,c("NLSI","CA")])



# find examplary points
pt3 <- SpatialPoints(sites[sites$NLSI > 0.20 &  sites$CA < 5000, ][1,c("X", "Y")]) #upper-left corner
pt2 <- SpatialPoints(sites[sites$NLSI > 0.15 &  sites$CA > 25000, ][1,c("X", "Y")]) #upper-right corner
pt1 <- SpatialPoints(sites[sites$NLSI < 0.05 &  sites$CA < 5000, ][1,c("X", "Y")]) #lower-left corner
pt4 <- SpatialPoints(sites[sites$NLSI < 0.07 &  sites$CA > 25000, ][1,c("X", "Y")]) #lower-right corner

pts <- rbind(pt3, pt2, pt1, pt4)
Buf <- gBuffer(pts, width=bufferWidth, byid = T)

bufferWidth <- 10000

par(mfrow=c(2,2), mar=c(2,2,2,2))
for(i in 1:length(Buf)){
  landBuf <- mask(CLC_SNH_NL, Buf[i])
  landBuf <- trim(landBuf, pad=2)
  plot(landBuf, legend = F, axes = F, box = F, col=c("gold1", "green4"))
}
