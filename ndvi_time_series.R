library(rts)
library(raster)
library(qdapRegex)
library(rgeos)
library(dplyr)
library(rio)
library(broom)

# load functions
source("functions.R")

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

## compute annual mean NDVI and annual sd
ndvi.annual.mean <- apply.yearly(ndvi.rts, mean)
ndvi.annual.sd <- apply.yearly(ndvi.rts, function(x)sd(x, na.rm=TRUE))

# optional: write or read data
data.ToReadWrite <- "ndvi.annual.mean" # set data name

dir.create(file.path(ndvi.Dir, "Processed_NDVI"), showWarnings = FALSE) # create dir

write.rts(get(data.ToReadWrite), paste0(ndvi.Dir,"/Processed_NDVI/", data.ToReadWrite), overwrite=T) # write rts
assign(data.ToReadWrite, read.rts(paste0(ndvi.Dir,"/Processed_NDVI/", data.ToReadWrite))) # read rts


##### extract by site (points) #####
# project points coordinates
coords <- spTransform(coords, CRS("+init=epsg:4326"))
sites[,2:3] <- coords@coords

### test different aggregation factors ###
# set aggregation factors
agr.fac <- 2:20

ndvi.points.sum <- c()
for(i in agr.fac){
  # aggregate to specified factor
  ndvi.annual.sd.agr <- aggregate(ndvi.annual.sd@raster, i)
  # extract data and append it
  ndvi.points.sum <- rbind.data.frame(ndvi.points.sum,
                                      cbind.data.frame(Site = sites[,1], agreg.fac = i,
                                                       sd.ndvi = apply(extract(ndvi.annual.sd.agr, sites[,2:3]), 1, 
                                                                       function(x)mean(x, na.rm = T))))
}

### extract data ###
# set aggregation factor
agr.fac <- 1

if(agr.fac > 1){
  ndvi.annual.sd.agr <- aggregate(ndvi.annual.sd@raster[[-nlayers(ndvi.annual.sd@raster)]], agr.fac)
  ndvi.annual.mean.agr <- aggregate(ndvi.annual.mean@raster[[-nlayers(ndvi.annual.mean@raster)]], agr.fac)
}else{
  ndvi.annual.sd.agr <- ndvi.annual.sd@raster[[-nlayers(ndvi.annual.sd@raster)]]
  ndvi.annual.mean.agr <- ndvi.annual.mean@raster[[-nlayers(ndvi.annual.mean@raster)]]
}

ndvi.points.sum <- cbind.data.frame(Site = sites[,1],
                                    stSD.ndvi = apply(extract(ndvi.annual.sd.agr, sites[,2:3]), 1, 
                                                      function(x)mean(x, na.rm = T)),
                                    ltSD.ndvi = apply(extract(ndvi.annual.mean.agr, sites[,2:3]), 1, 
                                                      function(x)sd(x, na.rm = T)),
                                    mean.ndvi = apply(extract(ndvi.annual.mean.agr, sites[,2:3]), 1, 
                                                      function(x)mean(x, na.rm = T)))


plot(ndvi.points.sum$mean.ndvi, ndvi.points.sum$stSD.ndvi) # show correlation btw. short-term and long-term measures
write.csv(ndvi.points.sum, paste0(ndvi.Dir,"/ndvi.points.sum.scv"), row.names = F) # write result of NDVI extraction
