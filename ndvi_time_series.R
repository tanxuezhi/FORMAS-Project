library(rts)
library(raster)
# library(lme4)
# library(lmerTest)
library(qdapRegex)
# library(visreg)
# library(cowplot)
library(rgeos)
# library(rgdal)
library(dplyr)


##### create 1km buffers #####
sites <- rio::import(file = "../Data/Butterflies - Netherlands/Sites.xlsx", which = 1L)
sites[,2:3] <- sites[,2:3]*1000
coords <- SpatialPoints(sites[,2:3], proj4string = CRS("+init=epsg:28992"))

buffer_NL <- gBuffer(coords, byid = T, width = 1000)
buffer_NL <- spTransform(buffer_NL, CRS("+init=epsg:4326"))

## write polygons to shp file
IDs <- sapply(slot(buffer_NL, "polygons"), function(x) slot(x, "ID"))
df <- data.frame(rep(0, length(IDs)), row.names=IDs)
SPDFxx <- SpatialPolygonsDataFrame(SPxx, df)
writeOGR(obj=buffer_NL, dsn=getwd(), layer="buffer_NL", driver="ESRI Shapefile")


##### import ndvi time serie #####

## download files from url list (generated from AppEEARS) 
url.List <- "../Data/NDVI/NDVI_NL/nl_250m_16days-download-list.txt"

## !!! doesn't work !!! ##
urls <- read.table(url.List)[1:209,]
for (url in urls) {
  download.file(url, destfile = basename(url))
}

##### process NDVI data #####

## directory containing NDVI files
ndvi.Dir <- "../Data/NDVI/NDVI_NL"

## extract dates from download list
time.ndvi <- read.table(url.List)
time.ndvi <- apply(time.ndvi, 1, as.character)
time.ndvi <- na.omit(unlist(ex_between(time.ndvi, "NDVI_doy", "_aid0001")))
time.ndvi <- as.Date(time.ndvi, format = "%Y%d%m")

## stack ndvi
ndvi <- stack(list.files(ndvi.Dir, pattern = ".tif", full.names = T))

## create raster time series
ndvi.rts <- rts(ndvi, time.ndvi)

## compute annual mean NDVI and annual CV
ndvi.annual.mean <- apply.yearly(ndvi.rts, mean)
ndvi.annual.sd <- apply.yearly(ndvi.rts, function(x)sd(x, na.rm=TRUE))

# optional: write or read data
data.ToReadWrite <- "ndvi.annual.mean" # set data name

dir.create(file.path(ndvi.Dir, "Processed_NDVI"), showWarnings = FALSE) # create dir

write.rts(get(data.ToReadWrite), paste0(ndvi.Dir,"/Processed_NDVI/", data.ToReadWrite), overwrite=T) # write rts
assign(data.ToReadWrite, read.rts(paste0(ndvi.Dir,"/Processed_NDVI/", data.ToReadWrite))) # read rts


##### extract by site (here buffers) #####

# extract annual mean values
ndvi.annual.mean.buffers <- c()
for (i in 1:nlayers(ndvi.annual.mean@raster)){
  ndvi.annual.mean.buffers.temp <- extract(ndvi.annual.mean[[i]], buffer_NL)
  ndvi.annual.mean.buffers <- rbind.data.frame(ndvi.annual.mean.buffers,
                                               cbind.data.frame(Site = names(buffer_NL), 
                                                                Year = ex_between(names(ndvi.annual.mean[[i]]), "X", ".")[[1]], 
                                                                mean.ndvi = unlist(lapply(ndvi.annual.mean.buffers.temp, function(x)mean(x, na.rm=T)))))
}

# extract annual sd values
ndvi.annual.sd.buffers <- c()
for (i in 1:nlayers(ndvi.annual.sd@raster)){
  ndvi.annual.sd.buffers.temp <- extract(ndvi.annual.sd[[i]], buffer_NL)
  ndvi.annual.sd.buffers <- rbind.data.frame(ndvi.annual.sd.buffers,
                                               cbind.data.frame(Site = names(buffer_NL), 
                                                                Year = ex_between(names(ndvi.annual.sd[[i]]), "X", ".")[[1]], 
                                                                sd.ndvi = unlist(lapply(ndvi.annual.sd.buffers.temp, function(x)sd(x, na.rm=T)))))
}

ndvi.annual.buffers <- cbind.data.frame(ndvi.annual.mean.buffers, sd.ndvi = ndvi.annual.sd.buffers[,3])
ndvi.annual.buffers$cv.ndvi <- ndvi.annual.buffers$sd.ndvi / ndvi.annual.buffers$mean.ndvi

# write data because of long computation time
save(ndvi.annual.buffers, file = paste0(ndvi.Dir,"/ndvi.annual.buffers.RData"))
load(paste0(ndvi.Dir,"/ndvi.annual.buffers.RData")) # reload data


# extract long-term and short-term CV
ndvi.buffers.sum <- ndvi.annual.buffers %>% 
  group_by(Site) %>% 
  summarise(LT.CV_ndvi = sd(mean.ndvi, na.rm=T)/mean(mean.ndvi, na.rm=T),
            mean_ndvi = mean(mean.ndvi, na.rm=T),
            CV.ST_ndvi = mean(cv.ndvi))

plot(ndvi.buffers.sum$LT.CV_ndvi, ndvi.buffers.sum$CV.ST_ndvi) # show correlation btw. short-term and long-term measures


##### analyses #####
lc_class <- rio::import(file = "//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Landcover/Corine_land-cover_2012_raster/Legend/clc_legend.xls", which = 1L)


cti_AB <- read.csv2("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Data/Butterflies - Netherlands/cti_site_year_abundance_2017-05-19.csv", dec = "."); cti_AB$Site <- as.factor(cti_AB$Site)
cti_AB <- merge(cti_AB, ndvi.site.final); cti_AB <- merge(cti_AB, sites); cti_AB$Landcover <- as.factor(cti_AB$Landcover)
cti_AB <- merge(cti_AB, lc_class, by.x = "Landcover", by.y = "CLC_CODE")

nb.year <- aggregate(cti ~ Site, cti_AB, length)

cti_AB <- cti_AB[cti_AB$Site %in% nb.year[nb.year[,2] > 10,1],]

cti_P <- read.csv2("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Data/Butterflies - Netherlands/cti_site_year_presence_2017-05-19.csv", dec = "."); cti_P$Site <- as.factor(cti_P$Site)
cti_P <- merge(cti_P, ndvi.site.final); cti_P <- merge(cti_P, sites); cti_P$Landcover <- as.factor(cti_P$Landcover)
cti_P <- merge(cti_P, lc_class, by.x = "Landcover", by.y = "CLC_CODE")

cti_P <- cti_P[cti_P$Site %in% nb.year[nb.year[,2] > 10,1],]


m_ltv_AB <- lmer(cti ~ Year * long.term.variability + (1|LABEL2/Site), data = cti_AB)
summary(m_ltv_AB)
pred.m_ltv_AB <- visreg(m_ltv_AB, xvar = "Year", ylab = "CTI", scale = "response", 
                        by = "long.term.variability", gg = T, breaks = 3, plot = F)

p.ltv_AB <- ggplot(data = pred.m_ltv_AB$fit, aes(x = Year, y = visregFit, color = as.factor(long.term.variability))) + 
  geom_line() + 
  scale_y_continuous(name = "CTI", limits=c(8.95,9.3)) +
  theme(legend.position = c(0.8,0.1)) + 
  ggtitle("Abundance-weighted") + 
  scale_color_manual(name = "NDVI variability", labels = c("Low", "Medium", "High"), values = c("royalblue1", "palegreen3", "tomato2")) 

m_ltv_PO <- lmer(cti ~ Year * long.term.variability + (1|LABEL2/Site), data = cti_P)
summary(m_ltv_PO)
pred.m_ltv_PO <- visreg(m_ltv_PO, xvar = "Year", ylab = "CTI", scale = "response", 
                        by = "long.term.variability", gg = T, breaks = 3, plot = F)

p.ltv_PO <- ggplot(data = pred.m_ltv_PO$fit, aes(x = Year, y = visregFit, color = as.factor(long.term.variability))) + 
  geom_line() + 
  scale_y_continuous(name = "CTI", limits=c(8.95,9.3)) +
  theme(legend.position = "none") + 
  ggtitle("Presence-only") + 
  scale_color_manual(name = "NDVI variability", labels = c("Low", "Medium", "High"), values = c("royalblue1", "palegreen3", "tomato2")) 

plot_grid(p.ltv_AB, p.ltv_PO)
