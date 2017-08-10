library(raster)
library(sp)
library(rgeos)
library(dplyr)
source("functions.R")

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

sites_SWE <- read.csv(file = "../Data/Birds - Sweden/Sites_SWE_ETRS89_landcover.csv")
sites_SWE <- SpatialPointsDataFrame(sites_SWE[,2:3], sites_SWE, proj4string = CRS("+init=epsg:3035"))
CLC_SNH_SWE <- crop(CLC_SNH, extent(sites_SWE) + 100000)
writeRaster(CLC_SNH_SWE, "../Landcover/Sweden/SNH_SWE.tif")

##### prepare Fragstat points #####
pts_Fragstat_FIN <- cbind.data.frame(ID=sites_FIN@data$Site, row = rowFromY(CLC_SNH_FIN, sites_FIN@coords[,"Y"]), col = colFromX(CLC_SNH_FIN, sites_FIN@coords[,"X"]))
pts_Fragstat_FIN <- data.frame(paste0("[",pts_Fragstat_FIN[,1],":",pts_Fragstat_FIN[,2],":",pts_Fragstat_FIN[,3],"]"))
colnames(pts_Fragstat_FIN) <- "FPT_TABLE"
write.table(pts_Fragstat_FIN, "../Connectivity/Fragmentation/FIN/FIN_sites.fpt", row.names = F, quote = F)

pts_Fragstat_NL <- cbind.data.frame(ID=sites_NL@data$Site, row = rowFromY(CLC_SNH_NL, sites_NL@coords[,"Y"]), col = colFromX(CLC_SNH_NL, sites_NL@coords[,"X"]))
pts_Fragstat_NL <- data.frame(paste0("[",pts_Fragstat_NL[,1],":",pts_Fragstat_NL[,2],":",pts_Fragstat_NL[,3],"]"))
colnames(pts_Fragstat_NL) <- "FPT_TABLE"
write.table(pts_Fragstat_NL, "../Connectivity/Fragmentation/NL/NL_sites.fpt", row.names = F, quote = F)

pts_Fragstat_SWE <- cbind.data.frame(ID=sites_SWE@data$Site, row = rowFromY(CLC_SNH_SWE, sites_SWE@coords[,"y"]), col = colFromX(CLC_SNH_SWE, sites_SWE@coords[,"x"]))

cbind.data.frame(pts_Fragstat_SWE[,1], 1:nrow(pts_Fragstat_SWE))

pts_Fragstat_SWE <- data.frame(paste0("[",1:nrow(pts_Fragstat_SWE),":",pts_Fragstat_SWE[,2],":",pts_Fragstat_SWE[,3],"]"))
colnames(pts_Fragstat_SWE) <- "FPT_TABLE"
write.table(pts_Fragstat_SWE, "../Connectivity/Fragmentation/SWE/SWE_sites.fpt", row.names = F, quote = F)


##### load results #####
folder <- "../Connectivity/Fragmentation/NL/"
dup.sites <- read.table("../Data/Butterflies - Netherlands/Duplicated_sites.txt", h = T)
frag <- extractFrag(folder, dup.sites, sites_NL)

# folder <- "../Connectivity/Fragmentation/FIN/"
# dup.sites <- read.table("../Data/Butterflies - Finland/Duplicated_sites.txt", h = T)
# frag <- extractFrag(folder, dup.sites, sites_FIN)

# folder <- "../Connectivity/Fragmentation/SWE/"
# dup.sites <- read.table("../Data/Birds - Sweden/Duplicated_sites.txt", h = T)
# frag <- extractFrag(folder, dup.sites, sites_SWE)

frag %>% group_by(Scale) %>% summarise(n = length(LID))

write.csv(frag, paste0(folder, "Frag_indices.csv"), row.names = F)


# corrplot(abs(cor(frag[,-c(1,2,36)])), method=c("color"),  
#          type="upper")
# plot(frag[,c("NLSI","CA","PROX_MN","CONNECT", "LSI")])
# plot(frag[,c("NLSI","CA")])


# find examplary points
pts <- read.csv("../Connectivity/Fragmentation/NL/Frag_indices.csv")
pts <- pts %>% filter(Scale == 5000)

sites_NL <- read.csv(file = "../Data/Butterflies - Netherlands/Sites_NL_ETRS89_landcover.csv")
sites_NL <- SpatialPointsDataFrame(sites_NL[,2:3], sites_NL, proj4string = CRS("+init=epsg:3035"))

plot(CLUMPY ~ PLAND, pts)
plot(NLSI ~ PLAND, pts)
plot(NLSI ~ CLUMPY, pts)


pt1 <- sites_NL[sites_NL$Site == pts[pts$CLUMPY < 0.80 &  pts$PLAND < 30, "LID"][1],]#lower-left corner
pt2 <- sites_NL[sites_NL$Site == pts[pts$CLUMPY < 0.80 &  pts$PLAND > 60, "LID"][1],] #lower-right corner
pt3 <- sites_NL[sites_NL$Site == pts[pts$CLUMPY > 0.90 &  pts$PLAND < 30, "LID"][1],] #upper-left corner
pt4 <- sites_NL[sites_NL$Site == pts[pts$CLUMPY > 0.90 &  pts$PLAND > 60, "LID"][1],] #upper-right corner

pts <- rbind(pt3, pt4, pt1, pt2)
Buf <- gBuffer(pts, width=5000, byid = T)

par(mfrow=c(2,2), mar=c(2,2,2,2))
for(i in 1:length(Buf)){
  landBuf <- mask(CLC_SNH_NL, Buf[i,])
  landBuf <- trim(landBuf, pad=2)
  plot(landBuf, legend = F, axes = F, box = F, col=c("gold1", "green4"))
}


# difference btw. countries
library(data.table)
library(ggplot2)
butterflies.data <- as.tbl(fread("../Data/cti_butterflies_data.csv"))
birds.data <- as.tbl(fread("../Data/cti_birds_data.csv"))

birds.data$country <- "SWE"

dat <- rbind(butterflies.data, birds.data)
dat <- subset(dat, dat$Scale == 5000)
dat <- dat %>% group_by (Site) %>% summarize_all(max)
# dat <- dat %>% filter(CLUMPY > .5)

ggplot(data = dat, aes(y = CLUMPY, x = PLAND, colour = country)) + geom_point() + theme_classic()

ggplot(data = dat, aes(y = PD, x = PLAND, colour = country)) + geom_point() + theme_classic() +
  facet_grid(~country)
