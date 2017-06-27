# load libraries
library(raster)
library(SDMTools)
library(rgeos)

#load sites
sites <- read.csv(file = "../Data/Butterflies - Netherlands/Sites_ETRS89_landcover.csv")
lc_class <- rio::import(file = "../Landcover/clc_legend.xls", which = 1L)
sites <- merge(sites, lc_class, by.x = "Landcover", by.y = "CLC_CODE")

#load landcover
CLC_NL <- raster("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Landcover/CLC_NL.tif")
# reclassify by second-level categories
reclassTable <- data.frame(from = lc_class[,1], to = as.numeric(as.factor(lc_class[,4])))
CLC_NL_reclas <- reclassify(CLC_NL, reclassTable)

#landuse stats at monitoring sites
table(sites$LABEL2)


#select land-use category
# Forest
LC <- 23:25

#reclass raster into habitat and non-habitat
CLC_NL_sel <- CLC_NL
CLC_NL_sel[CLC_NL_sel %in% LC] <- 1
CLC_NL_sel[!CLC_NL_sel %in% 1] <- 0

#create buffers
pts <- sites[!sites$LABEL1 %in% "Artificial surfaces",]

landFrag <- c()
for(i in c(1:nrow(pts))){
  landPt <- reclassTable[reclassTable$from %in% pts[i,"GRID_CODE"],2]
  
  buf <- gBuffer(SpatialPoints(pts[i,3:4]), width = 3000)
  cr <- crop(CLC_NL_reclas, extent(buf), snap="out")
  fr <- rasterize(buf, cr)
  lr <- mask(x=cr, mask=fr)
  
  stat <- ClassStat(lr)
  stat <- stat[stat$class %in% landPt,]

  landFrag <- rbind.data.frame(landFrag, cbind.data.frame(Site = pts[i,2], stat))
}

write.csv(landFrag, "../Data/Butterflies - Netherlands/landFrag.csv", row.names = F)

