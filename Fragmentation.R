library(corrplot)


###load landcover
# reclassified by semi-natural habitat
CLC_NL_SNH <- raster("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Landcover/Corine_land-cover_2012_raster/CLC_2012_100m_reclass_snh_NL.tif")

###load sites
sites <- read.csv(file = "../Data/Butterflies - Netherlands/Sites_ETRS89_landcover.csv")

##### prepare Fragstat points #####
pts_Fragstat <- cbind.data.frame(ID=sites[,1], row = rowFromY(CLC_NL_SNH, sites[,3]), col = colFromX(CLC_NL_SNH, sites[,2]))
data.frame(paste0("[",pts_Fragstat[,1],":",pts_Fragstat[,2],":",pts_Fragstat[,3],"]"))

##### load results #####
frag <- read.table("../Connectivity/Fragmentation/NL/results.class", h=T,na.strings="N/A")

corrplot(abs(cor(frag[,-c(1,2,36)])), method=c("color"),  
         type="upper")
plot(frag[,c("NLSI","CA","PROX_MN")])

frag[,1] <- as.numeric(gsub("point_", "", frag[,1]))
frag <- frag[,c("LID","NLSI","CA","PROX_MN")]

write.csv(frag, "../Connectivity/Fragmentation/NL/Frag_indices_NL.csv", row.names = F)
