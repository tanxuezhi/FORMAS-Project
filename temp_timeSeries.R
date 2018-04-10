library(raster)
library(rts)
library(stringr)
library(ggplot2)
library(data.table)


library(fmi)
apiKey <- "2da0d2bb-73b0-4b58-9200-8f97cf5645c1"
request <- FMIWFSRequest$new(apiKey = apiKey)

request$setParameters(request = "getFeature",
                      storedquery_id = "fmi::observations::weather::monthly::grid",
                      starttime = "2013-01-01T00:00:00Z",
                      endtime = "2016-01-01T00:00:00Z",
                      bbox = "19.09,59.3,31.59,70.13")
client <- FMIWFSClient$new(request = request)
response <- client$getRaster(parameters = list(splitListFields = TRUE))

##### Finland ####

r <- brick("C:/Local Folder (c)/monthly_19610101T000000_19900101T000000_20161201T000000_latlon.grb2")
names.layers <- paste0(rep(c(1990:2016), each = 12))

library(rNOMADS)
r2 <- ReadGrib("C:/Local Folder (c)/monthly_19610101T000000_19900101T000000_20161201T000000_latlon.grb2")


##### European scale #####
r <- brick("//storage.slu.se/Home$/yofo0001/Desktop/tg_0.25deg_reg_v15.0.nc")
r <- subset(r, names(r)[grep(c("\\.03\\.|\\.04\\.|\\.05\\.|\\.06\\.|\\.07\\.|\\.08\\."), names(r))])
r <- subset(r, names(r)[as.numeric(str_sub(names(r), 2,5)) > 1988])

rasterTS <- rts(r, as.Date(gsub("X", "", names(r)), format = "%Y.%m.%d"))
rasterTSMonth <- apply.yearly(rasterTS, mean)
write.rts(rasterTSMonth, "../Temperature/rasterTSMonth", overwrite = T)


####

rasterTSMonth <- read.rts("../Temperature/rasterTSMonth")

birds.data <- as.tbl(fread("../Data/cti_birds_data.csv")) %>% dplyr:::filter(type == "Presence", Scale == 3000)

butterflies.data <- as.tbl(fread("../Data/cti_butterflies_data.csv"))
butterflies.data <- butterflies.data %>% dplyr:::filter(type == "Presence", Scale == 3000)
butterflies.data <- 
  left_join(
    butterflies.data %>%
      group_by(coords = paste(X,Y)) %>% summarise(Site = first(Site)) %>%
      ungroup(),
    butterflies.data %>% mutate(coords = paste(X,Y)) %>% dplyr:::select(-2),
    by = c("coords" = "coords")
  ) %>% dplyr:::select(-1) %>%
  group_by(Site, Year) %>% summarise_all(function(x){ifelse(is.numeric(x),mean(x, na.rm = TRUE), x)})


cti.data <- butterflies.data
trend <- c()
n = 0
for(i in 1:length(unique(cti.data$Site))){
  n = n+1
  
  site <- cti.data %>% dplyr:::filter(Site == unique(cti.data$Site)[i])
  Coords <- SpatialPoints(site[1,c("X", "Y")], proj4string = CRS("+init=epsg:3035"))
  Coords <- spTransform(Coords, CRS("+init=epsg:4326"))
  time <- 1990:2016
  time <- paste0(time, "0831")
  
  temp_sites <- as.data.frame(extract(rasterTSMonth, Coords, time))
  temp_sites$Year <- as.numeric(substr(rownames(data.frame(temp_sites)), 1,4))
  colnames(temp_sites)[1] <- "Temp" ; rownames(temp_sites) <- NULL
  
  if(any(is.na(temp_sites[1]))){
    rasterTSMonth.ag <- rasterTSMonth
    rasterTSMonth.ag@raster <- aggregate(rasterTSMonth.ag@raster, 2)
    temp_sites <- as.data.frame(extract(rasterTSMonth.ag, Coords, time))
    temp_sites$Year <- as.numeric(substr(rownames(data.frame(temp_sites)), 1,4))
    colnames(temp_sites)[1] <- "Temp" ; rownames(temp_sites) <- NULL
    
  }
  
  if(any(is.na(temp_sites[1]))){
    rasterTSMonth.ag <- rasterTSMonth
    rasterTSMonth.ag@raster <- aggregate(rasterTSMonth.ag@raster, 4)
    temp_sites <- as.data.frame(extract(rasterTSMonth.ag, Coords, time))
    temp_sites$Year <- as.numeric(substr(rownames(data.frame(temp_sites)), 1,4))
    colnames(temp_sites)[1] <- "Temp" ; rownames(temp_sites) <- NULL
    
  }
  
  if(any(is.na(temp_sites[1]))){
    rasterTSMonth.ag <- rasterTSMonth
    rasterTSMonth.ag@raster <- aggregate(rasterTSMonth.ag@raster, 8)
    temp_sites <- as.data.frame(extract(rasterTSMonth.ag, Coords, time))
    temp_sites$Year <- as.numeric(substr(rownames(data.frame(temp_sites)), 1,4))
    colnames(temp_sites)[1] <- "Temp" ; rownames(temp_sites) <- NULL
    
  }
  
  trend.cti <- lm(cti ~ Year, site)
  trend.t <- lm(Temp ~ Year, temp_sites)
  
  trend.temp <- cbind.data.frame(Site = unique(cti.data$Site)[i], 
                                 nYears = nrow(site),
                                 trend.temp = trend.t$coefficients[2],
                                 trend.cti = trend.cti$coefficients[2])
  trend.temp$trend.dif <- trend.temp$trend.cti - trend.temp$trend.t
  trend.temp$trend.ratio <- trend.temp$trend.cti / trend.temp$trend.t
  
  trend <- rbind.data.frame(trend, trend.temp)
  rownames(trend) <- NULL
  
  print(paste0("Site ",unique(cti.data$Site)[i], " : ",  
               round((n/length(unique(cti.data$Site)))*100, digits = 4), "%"))
  
  write_csv(trend.temp, "../cti.temp.trends.csv", append = T, col_names = F)
}

################
## write data ##
################

#birds
cti.trends <- fread("../Data/cti.temp.trends_birds.csv")
cti.trends <- as.tbl(cti.trends)

birds.data <- as.tbl(fread("../Data/cti_birds_data.csv"))
birds.data <- birds.data %>% dplyr:::filter(Scale %in% c(3000,30000), type == "Presence")
birds.data <- birds.data %>% group_by(Site, Scale) %>% dplyr:::filter(row_number(type) == 1)
birds.data <- birds.data[,c(3,5,6,8,10,11,17,18)]

cti.trends.data <- left_join(cti.trends, birds.data, by = c("Site"))

write_csv(cti.trends.data, "../Data/cti.trends_birds_data.csv")


#butterflies
cti.trends <- read_csv("../Data/cti.temp.trends_butterflies.csv")

butterflies.data <- as.tbl(fread("../Data/cti_butterflies_data.csv"))
butterflies.data <- butterflies.data %>% dplyr:::filter(Scale %in% c(3000,30000), type == "Presence")
butterflies.data <- butterflies.data %>% group_by(Site, Scale) %>% dplyr:::filter(row_number(type) == 1)
butterflies.data <- butterflies.data[,c(2,5,6,8,10,11,17,18,19)]

cti.trends.data <- left_join(cti.trends, butterflies.data, by = c("Site"))

write_csv(cti.trends.data, "../Data/cti.trends_butterflies_data.csv")
