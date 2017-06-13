library(rts)
library(raster)
# library(lme4)
# library(lmerTest)
library(qdapRegex)
# library(visreg)
# library(cowplot)
library(rgeos)
# library(rgdal)

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
ndvi.annual.cv <- apply.yearly(ndvi.rts, function(x)sd(x, na.rm=TRUE) / mean(x, na.rm=TRUE))


annual.cv.ndvi <- calc(ndvi.annual@raster, function(x)sd(x, na.rm=TRUE)) / calc(ndvi.annual@raster, function(x)mean(x, na.rm=TRUE))

# optional: write or read data
dir.create(file.path(ndvi.Dir, "Processed_NDVI"), showWarnings = FALSE) # create dir
data.ToReadWrite <- "ndvi.annual.mean"

write.rts(get(data.ToReadWrite), paste0(ndvi.Dir,"/Processed_NDVI/", data.ToReadWrite), overwrite=T) # write rts
assign(data.ToReadWrite, read.rts(paste0(ndvi.Dir,"/Processed_NDVI/", data.ToReadWrite))) # read rts



##### extract by site #####

coords <- spTransform(coords, CRS("+init=epsg:4326"))
sites[,2:3] <- coords@coords


ndvi.ts.pts <- extract(ndvi.rts, coords)

ndvi.site = list()
for (i in 1:ncol(ndvi.ts.pts)){
  
  temp <- ndvi.ts.pts[,i]
  cv.temp <- apply.yearly(temp, function(x)sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE))
  mean.temp <- apply.yearly(temp, function(x)mean(x, na.rm =T))
  
  temp.res <- cbind.data.frame(Site = sites[i,1], 
                               year = substr(rownames(data.frame(mean.temp)), 1,4), 
                               mean.ndvi = mean.temp, 
                               cv.ndvi = cv.temp)
  row.names(temp.res) <- NULL
  
  ndvi.site[[i]] <- temp.res
}
ndvi.site = do.call(rbind, ndvi.site)

ndvi.site.final <- aggregate(mean.ndvi ~ Site, ndvi.site, function(x)sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE))
ndvi.site.final$yearly.variability <- aggregate(cv.ndvi ~ Site, ndvi.site, function(x)mean(x, na.rm =T))[,2]
colnames(ndvi.site.final) <- c("Site", "long.term.variability", "yearly.variability")
ndvi.site.final$mean.ndvi <- na.omit(apply(ndvi.ts.pts,2,function(x)mean(x, na.rm =T)))



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


library(effects)
effect("long.term.variability", m_ltv_AB)

#### test with yearly variability ####
m_yv <- lmer(cti ~ Year * yearly.variability + (1|LABEL2/Site), data = cti_AB)
summary(m_yv)
visreg(m_yv, xvar = "Year", ylab = "CTI", scale = "response", ylim = c(8.9,9.3), 
       by = "yearly.variability", gg = T, breaks = 3)



##### test correlations #####
library(HH)

sit <- merge(sites, ndvi.site.final); sit <- merge(sit, lc_class, by.x = "Landcover", by.y = "CLC_CODE")


vif(sit[,c("long.term.variability", "Landcover", "yearly.variability")])

