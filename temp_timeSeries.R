library(raster)
library(rts)
library(stringr)
library(ggplot2)
library(data.table)
library(dplyr)

r <- brick("//storage.slu.se/Home$/yofo0001/Desktop/tg_0.25deg_reg_v15.0.nc")
r <- subset(r, names(r)[grep(c("\\.03\\.|\\.04\\.|\\.05\\.|\\.06\\.|\\.07\\.|\\.08\\."), names(r))])
r <- subset(r, names(r)[as.numeric(str_sub(names(r), 2,5)) > 1988])

rasterTS <- rts(r, as.Date(gsub("X", "", names(r)), format = "%Y.%m.%d"))
rasterTSMonth <- apply.yearly(rasterTS, mean)
write.rts(rasterTSMonth, "../Temperature/rasterTSMonth", overwrite = T)

rasterTSMonth <- read.rts("../Temperature/rasterTSMonth")


birds.data <- as.tbl(fread("../Data/cti_birds_data.csv"))
birds.data <- birds.data %>% filter(Scale == 5000, type == "Abundance")

butterflies.data <- as.tbl(fread("../Data/cti_butterflies_data.csv"))
butterflies.data <- butterflies.data %>% filter(Scale == 5000, type == "Abundance")


trend <- c()
for(i in 1:length(unique(butterflies.data$Site))){
  site <- butterflies.data %>% filter(Site == unique(butterflies.data$Site)[i])
  Coords <- SpatialPoints(site[1,c("X", "Y")], proj4string = CRS("+init=epsg:3035"))
  Coords <- spTransform(Coords, CRS("+init=epsg:4326"))
  time <- range(site$Year)[1]:range(site$Year)[2]
  time <- paste0(time, "0831")
  temp_sites <- extract(rasterTSMonth, Coords, time)
  
  trend.cti <- lm(cti ~ c(1:nrow(site)), site)
  
  
  if(!is.na(temp_sites[1])){
    trend.temp <- cbind.data.frame(Site = unique(butterflies.data$Site)[i], 
                                   trend.t = lm(temp_sites ~ c(1:length(temp_sites)))$coefficients[2],
                                   trend.cti = lm(cti ~ c(1:nrow(site)), site)$coefficients[2])
    trend.temp$debt <- trend.temp$trend.cti - trend.temp$trend.t
  }else{
    trend.temp <- cbind.data.frame(Site = unique(butterflies.data$Site)[i], trend.t = NA,
                                   trend.cti = lm(cti ~ c(1:nrow(site)), site)$coefficients[2],
                                   debt = NA)
  }
  trend <- rbind.data.frame(trend, trend.temp)
  rownames(trend) <- NULL
}

write.csv(trend, "../Temperature/trend_FIN_NL.csv", row.names = F)

trend <- read.csv("../Temperature/trend_SWE.csv")


temp <- c()
for(i in 1:length(unique(birds.data$Site))){
  site <- birds.data %>% filter(Site == unique(birds.data$Site)[i])
  Coords <- SpatialPoints(site[1,c("X", "Y")], proj4string = CRS("+init=epsg:3035"))
  Coords <- spTransform(Coords, CRS("+init=epsg:4326"))
  temp_sites <- extract(rasterTSMonth, Coords, paste0(site$Year, "0831"))
  
  temp <- rbind.data.frame(temp, cbind.data.frame(Site = unique(birds.data$Site)[i],
                                                  cti = site$cti, 
                                                  X = site$X,
                                                  Y = site$Y,
                                                  temp = temp_sites))
  rownames(temp) <- NULL
}

library(lmerTest)
m <- lmer(cti ~ temp + X + Y + (1 + temp|Site), data = temp)
summary(m)
library(visreg)
par(mfrow=c(2,2))
visreg(m, xvar = c("temp"))

library(broom)
coef.slope <- na.omit(temp) %>% group_by(Site) %>% do(tidy(lm(cti ~ temp, data = .)))
coef.slope <- coef.slope %>% filter(term == "temp")


m2 <- lm(temp ~ Y, temp)
summary(m2)
visreg(m2)
