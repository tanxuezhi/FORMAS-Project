library(corrplot)
library(rgeos)

###load landcover
# reclassified by semi-natural habitat
CLC_SNH_NL <- raster("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Landcover/Netherlands/CLC_2012_100m_reclass_snh_NL.tif")

###load sites
sites <- read.csv(file = "../Data/Butterflies - Netherlands/Sites_ETRS89_landcover.csv")

##### prepare Fragstat points #####
pts_Fragstat <- cbind.data.frame(ID=sites[,1], row = rowFromY(CLC_NL_SNH, sites[,3]), col = colFromX(CLC_NL_SNH, sites[,2]))
data.frame(paste0("[",pts_Fragstat[,1],":",pts_Fragstat[,2],":",pts_Fragstat[,3],"]"))

##### load results #####
frag <- read.table("../Connectivity/Fragmentation/NL/results_10000radius_3000connect.class", h=T,na.strings="N/A")

corrplot(abs(cor(frag[,-c(1,2,36)])), method=c("color"),  
         type="upper")
plot(frag[,c("NLSI","CA","PROX_MN","CONNECT", "LSI")])
plot(frag[,c("NLSI","CA")])


frag[,1] <- as.numeric(gsub("point_", "", frag[,1]))
frag <- frag[,c("LID","NLSI","CA","PROX_MN", "LSI", "ED")]

write.csv(frag, "../Connectivity/Fragmentation/NL/Frag_indices_NL.csv", row.names = F)



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
