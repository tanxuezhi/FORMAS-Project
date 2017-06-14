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
ndvi.annual.mean.buffers <- extractRTS(ndvi.annual.mean, buffer_NL, function(x)mean(x, na.rm =T))
ndvi.annual.sd.buffers <- extractRTS(ndvi.annual.mean, buffer_NL, function(x)sd(x, na.rm =T))

ndvi.annual.buffers <- cbind.data.frame(ndvi.annual.mean.buffers[,1:2], 
                                        mean.ndvi = ndvi.annual.mean.buffers[,3], 
                                        sd.ndvi = ndvi.annual.sd.buffers[,3])
ndvi.annual.buffers$cv.ndvi <- ndvi.annual.buffers$sd.ndvi / ndvi.annual.buffers$mean.ndvi
ndvi.annual.buffers <- filter(ndvi.annual.buffers, Year < 2017)

# detrend mean ndvi
trends <- ndvi.annual.buffers %>% 
  na.omit() %>%
  group_by(Site) %>% 
  do(trends = lm(mean.ndvi ~ Year, .))

detrended.ndvi <- select(trends %>% augment(trends), Site, Year, .std.resid)

ndvi.annual.buffers <- merge(ndvi.annual.buffers, detrended.ndvi)
names(ndvi.annual.buffers)[6] <- "detrended.ndvi"

ggplot(ndvi.annual.buffers, aes(x = Year, y = detrended.ndvi, color = Site)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F) +
  theme(legend.position = "none")

# write data because of long computation time
save(ndvi.annual.buffers, file = paste0(ndvi.Dir,"/ndvi.annual.buffers.RData"))
load(paste0(ndvi.Dir,"/ndvi.annual.buffers.RData")) # reload data


# extract long-term and short-term CV
ndvi.buffers.sum <- ndvi.annual.buffers %>% 
  group_by(Site) %>% 
  summarise(sd.LT_ndvi = sd(detrended.ndvi),
            mean_ndvi = mean(mean.ndvi, na.rm=T),
            CV.ST_ndvi = mean(cv.ndvi))

plot(ndvi.buffers.sum$sd.LT_ndvi, ndvi.buffers.sum$CV.ST_ndvi) # show correlation btw. short-term and long-term measures

write.csv(ndvi.buffers.sum, paste0(ndvi.Dir,"/ndvi.buffers.sum.scv"), row.names = F) # write result of NDVI extraction
