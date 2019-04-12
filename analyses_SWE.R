library(tidyverse)
library(lmerTest)
library(spatialEco)
library(BBmisc)
library(raster)
library(gstat)
library(nlme)

### load data and compute cti
data <- read.csv("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/Temporal_climatic_debt/SWE_PA_data.csv") %>% 
  filter(latin != "Branta leucopsis", yr > 1998) %>% as.tbl

frag_data <- read_csv("../Connectivity - Fragmentation/Fragmentation/frag_SWE.csv")
sti <- read_csv("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Data/Birds Atlas/EBCC1/EBCC1_sti.csv")

data <- data %>% rename(species = latin) %>% left_join(sti) %>%
  filter(abundance > 0) %>% group_by(karta, yr) %>%
  summarise(CTI_ab = sum(sti * abundance) / sum(abundance), CTI_pres = mean(sti)) %>%
  left_join(frag_data)

sites.SWE <- rio::import("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/FORMAS project/Data/Birds - Sweden/Koordinater standardrutternas samtliga punkter.xlsx")
sites.SWE <- sites.SWE %>% group_by(karta) %>% summarise(x = mean(rt90_o), y = mean(rt90_n)) 

### load temperature
tmp <- stack("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/CRU_data\\breeding\\tmp.grd")
names(tmp) <- 1901:2017
crs(tmp) <- crs("+init=epsg:3021")

### sweden mask
SWE <- raster::getData('GADM', country='SWE', level=0)
SWE <- spTransform(SWE, CRSobj = crs(tmp))

### compute trends
msk_sites <- rasterFromXYZ(sites.SWE[,-1] %>% 
                             complete(x,y),crs = crs(tmp))

lm_fun = function(x){
  if(length(x[!is.na(x)]) < 2){
    NA
  } else {
    df <- cbind.data.frame(y = as.vector(x), time = 1:length(x))
    m = gls(y ~ time, data = na.omit(df), correlation = corAR1(form = ~ time))
    return(summary(m)$coefficients[2])
  }
}


meanDiff_fun = function(x){
  if(length(x[!is.na(x)]) < 2){
    NA
  } else {
    df <- cbind.data.frame(y = as.vector(x), time = 1:length(x))
    df <- na.omit(df)
    df <- df %>% mutate(diff.y = y - lag(y), diff.time = time - lag(time), diff.y = diff.y / diff.time)
    return(mean(df$diff.y, na.rm = T))
  }
}


## cti
cti_raster <- vector("list", length(unique(data$yr)))
for(i in 1:length(unique(data$yr))){
  cti_raster[[i]] <- merge(rasterFromXYZ(data[data$yr == unique(data$yr)[i], c(3:6)] %>% 
                                           complete(x,y),crs = crs(tmp))[[2]], msk_sites)
}
cti_raster <- stack(cti_raster)
names(cti_raster) <- unique(data$yr)

ctiTrend <- calc(cti_raster, lm_fun)
plot(ctiTrend)

# interpolate
trend.df <- na.omit(as.data.frame(ctiTrend, xy=T))
coordinates(trend.df) =~ x+y
ctiTrend.idw <- autoKrige(layer ~ 1, trend.df)
ctiTrend.idw <- raster(ctiTrend.idw$krige_output)
proj4string(ctiTrend.idw) <- proj4string(tmp)
ctiTrend.idw <- mask(ctiTrend.idw, SWE)

## temp
tmp_sub <- tmp[[98:117]]
tmpTrend <- calc(tmp_sub, lm_fun)
plot(tmpTrend)

# interpolate
tmpTrend.df <- na.omit(as.data.frame(tmpTrend, xy=T))
coordinates(tmpTrend.df) =~ x+y
tmpTrend.idw <- autoKrige(layer ~ 1, tmpTrend.df)
tmpTrend.idw <- raster(tmpTrend.idw$krige_output)
proj4string(tmpTrend.idw) <- proj4string(tmp)
tmpTrend.idw <- mask(tmpTrend.idw, SWE)

# plot
# raw
par(mfrow=c(1,2))
plot(tmpTrend)
plot(ctiTrend)
par(mfrow=c(1,1))

# interpolated
par(mfrow=c(1,2))
plot(tmpTrend.idw)
plot(ctiTrend.idw)
par(mfrow=c(1,1))
 
# neigh = 8
# 
# library(doFuture)
# registerDoFuture()
# plan(multiprocess, workers = 3)
# 
# trends_sliding_window <- foreach(i = unique(data$karta), .combine = rbind.data.frame) %dopar% {
#   
#   focal_site <- data %>% filter(karta == i)
#   
#   regional_sites <- sites.SWE[pointDistance(unique(focal_site[, c("x","y")]), sites.SWE[, c("x","y")], lonlat = F) < user_dist,]
#   
#   regional_data <- data %>% filter(karta %in% regional_sites$karta)
#   
#   regional_cti_trend <- lme(CTI_ab ~ yr, random = list(karta =  ~ 1),
#                             correlation = corExp(form = ~ jitter(x) + jitter(y)), data = regional_data)
#   
#   regional_cti_trend <- regional_cti_trend$coefficients$fixed[[2]]
#   
#   regional_temp <- cbind.data.frame(regional_sites, raster::extract(tmp, regional_sites[,-1]))
#   
#   regional_temp <- regional_temp %>% gather(-1,-2,-3, key = "yr", value = "temp") %>% 
#     mutate(yr = as.numeric(gsub("layer.", "", yr)) + 1900)
#   
#   
#   regional_temp_trend <- lme(temp ~ yr, random = list(karta =  ~ 1),
#                              correlation = corExp(form = ~ jitter(x) + jitter(y)), data = regional_temp)
#   
#   regional_temp_trend <- regional_temp_trend$coefficients$fixed[[2]]
#   
#   
#   cbind.data.frame(unique(focal_site[,c(1,5:8)]), 
#                    regional_temp_trend, regional_cti_trend,
#                    nSites = length(unique(regional_sites$karta)))
# }
# 
# 
# plot(regional_cti_trend ~ regional_temp_trend, data = trends_sliding_window)
# 
# plot(regional_cti_trend ~ nSites, data = trends_sliding_window)
# 
# 
# plot(rasterFromXYZ(trends_sliding_window[,-1] %>% complete(x,y)))
# 
# trends_sliding_window <- trends_sliding_window %>% filter(!nSites < 6) %>%
#   mutate(CTA = (max((scale(regional_cti_trend) - scale(regional_temp_trend))) - 
#                   (scale(regional_cti_trend) - scale(regional_temp_trend))) * 
#            (scale(regional_cti_trend) + scale(regional_temp_trend)))
# 
# 
# trends_sliding_window$CTA <- normalize(as.vector(trends_sliding_window$CTA), method = "range")
# 
# 
# m.CTA <- gam(CTA ~ CLUMPY * PLAND + te(x, y),
#              data = trends_sliding_window)
# summary(m.CTA)



# analysis

# cti matching

trend <- cbind.data.frame(sites.SWE, raster::extract(stack(tmpTrend, ctiTrend), sites.SWE[,-1])) %>% as.tbl
names(trend)[4:5] <- c("tmpTrend", "ctiTrend")

trend <- data %>% group_by(karta) %>%
  summarise(nYr = length(unique(yr))) %>% mutate(tmp = raster::extract(mean(tmp), unique(data[,c("x","y")]))) %>% 
  left_join(trend)

trend <- trend %>% mutate(CTA = (max((scale(ctiTrend) - scale(tmpTrend))) - 
                                   (scale(ctiTrend) - scale(tmpTrend))) * 
                            (scale(ctiTrend) + scale(tmpTrend)))
trend$CTA <- normalize(as.vector(trend$CTA), method = "range")

trend <- trend %>% left_join(frag_data)

but.pts <- SpatialPoints(trend[,c("x", "y")])
but.r <- raster(ext = extent(but.pts)*1.5, resolution = 300000)
values(but.r) <- c(1:ncell(but.r))
gridCell <- cbind.data.frame(trend[,c(1,4,5)], gridCell = extract(but.r, but.pts))

trend <- trend %>% left_join(gridCell[,c(1,4)], by = "karta")

plot(rasterFromXYZ(trend[,-1] %>% complete(x,y)))
plot(nYr ~ ctiTrend, data = trend); abline(v = 0)
plot(CTA ~ tmpTrend, data = trend)



library(gamm4)
m.trend <- lmer(CTA ~ CLUMPY * PLAND + tmp + x*y + (1|gridCell),
               weight = nYr,
              data = trend)
summary(m.trend)

library(visreg)
visreg(m.trend, xvar = "CLUMPY", by = "PLAND")

data2 <- data %>% left_join(data %>% group_by(karta) %>% summarise(x = unique(x), y = unique(y)) %>% 
                              mutate(tmp = raster::extract(mean(tmp), unique(data[,c("x","y")])))) %>% left_join(gridCell)

m.trend <- gamm(CTI_pres ~ yr * CLUMPY * PLAND + tmp + te(x,y), 
                random = list(karta = ~ 1), 
                correlation = corAR1(form = ~ yr | karta),
               data = data2)
summary(m.trend$gam)



plot(ctiTrend ~ PLAND, data = trend)
plot(ctiTrend ~ CLUMPY, data = trend)
plot(ctiTrend ~ tmp, data = trend)


fitDist(fullData$CTI_ab)

m.ab <- gamlss(CTI_ab ~ CLUMPY * PLAND * yr + x*y + random(karta), family = c("SHASHo", "Sinh-Arcsinh"),
               data = fullData %>% ungroup %>% mutate(karta = as.factor(karta)))
plot(m.ab)
summary(m.ab)

m.ab <- lmer(log(CTI_ab) ~ CLUMPY * PLAND * yr + x*y + (yr|karta), 
             data = fullData)
car::qqPlot(resid(m.ab))
hist(resid(m.ab))
summary(m.ab)

m.pres <- lmer(CTI_pres ~ CLUMPY * PLAND *  yr + x*y + (yr|karta), 
               data = fullData)
summary(m.pres)
