#load landcover
CLC_NL <- raster("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Landcover/CLC_NL.tif")

#select land-use category
CLC_NL_forest <- CLC_NL
CLC_NL_forest[CLC_NL_forest %in% 23:25] <- 1
CLC_NL_forest[!CLC_NL_forest %in% 1] <- 0

# aggregate raster
CLC_NL_forest <- aggregate(CLC_NL_forest, 4, max)
