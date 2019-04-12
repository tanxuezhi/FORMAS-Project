library(tidyverse)
library(raster)
library(data.table)
library(rgeos)
library(emmeans)
source("functions.R")

#####################
### Illustrations ###
#####################

##############
## data set ##
##############

butterflies.data.presence <- as.tbl(fread("../Data/cti_butterflies_data.csv")) %>% dplyr::filter(type == "Presence")
butterflies.data.presence %>% dplyr::filter(Scale == 50000) %>% group_by(country) %>% summarise(n = length(unique(Site)))
butterflies.data.presence %>% dplyr::filter(Scale == 50000) %>% group_by(country) %>% summarise(first = min(Year), last = max(Year))

butterflies.data.presence %>% dplyr::filter(Scale == 50000) %>% group_by(country, Site) %>% summarise(n = length(unique(Year))) %>%
  group_by(country) %>%
  summarise(mean.noYear = mean(n), sd.noYear = sd(n))

butterflies.data.presence %>% dplyr::filter(Scale == 50000) %>% group_by(country, Year) %>% summarise(n = length(unique(Site))) %>%
  group_by(country) %>%
  summarise(mean.noSite = mean(n), sd.noSite = sd(n))


data <- butterflies.data.presence %>% dplyr::filter(Scale == 50000) %>% 
  group_by(Site) %>% summarise(nYear = n(), X = unique(X), Y = unique(Y))

ggplot(data, aes(x = X, y = Y, color = nYear)) + geom_point(alpha = .5) + 
  scale_color_gradientn(name = "No. years \n", 
                       breaks =  c(1,max(data$nYear)),
                       limits = c(1,max(data$nYear)),
                       colours = c("yellow","red")) +
  scale_x_continuous("Longitude") + scale_y_continuous("Latitude")


###############
## landscape ##
###############

butterflies.data.presence %>% dplyr::filter(Scale == 50000) %>% group_by(country, Year) %>% Hmisc::describe()
  
ggplot(butterflies.data.presence %>% dplyr:::filter(PLAND < .95 & PLAND > .05) %>% group_by(country, Scale, Habitat) %>%
         distinct(), aes(y = CLUMPY, x = PLAND, color =country)) + 
  geom_point(alpha = .3) + facet_wrap( ~ Scale, labeller = label_wrap) + 
  scale_color_manual("Country", values = c("#00008B", "#BDB76B"), labels = c("Finland", "Netherlands")) + 
  scale_y_continuous("Habitat clumpiness") + scale_x_continuous("% Semi-natural habitat")

SNH_open <- raster("../Landcover/SNH/SNH_open.tif")
SNH_forest <- raster("../Landcover/SNH/SNH_forest.tif")

butterflies.data.presence %>% dplyr::filter(Scale == 20000, CLUMPY < .8, PLAND > .7)

# high clump, high snh
pt1 <- butterflies.data.presence %>% dplyr::filter(Site == "1231_NL") %>% dplyr::select(X, Y) %>% summarise_all(unique)
buf1 <- gBuffer(spgeom = SpatialPoints(pt1), width = 20000)
land1 <- raster::mask(crop(SNH_open, buf1), buf1)

# high clump, low snh
pt2 <- butterflies.data.presence %>% dplyr::filter(Site == "1743_NL") %>% dplyr::select(X, Y) %>% summarise_all(unique)
buf2 <- gBuffer(spgeom = SpatialPoints(pt2), width = 20000)
land2 <- raster::mask(crop(SNH_open, buf2), buf2)

# low clump, high snh
pt3 <- butterflies.data.presence %>% dplyr::filter(Site == "72_FIN") %>% dplyr::select(X, Y) %>% summarise_all(unique)
buf3 <- gBuffer(spgeom = SpatialPoints(pt3), width = 20000)
land3 <- raster::mask(crop(SNH_forest, buf3), buf3)

# low clump, low snh
pt4 <- butterflies.data.presence %>% dplyr::filter(Site == "77_FIN") %>% dplyr::select(X, Y) %>% summarise_all(unique)
buf4 <- gBuffer(spgeom = SpatialPoints(pt4), width = 20000)
land4 <- raster::mask(crop(SNH_open, buf4), buf4)

butterflies.data.presence %>% dplyr::filter(Site %in% c("1231_NL", "1743_NL", "72_FIN", "77_FIN"), Scale == 20000, Year == 2010)

# plot
par(mfrow=c(2,2), mar = c(0.1,0.1,0.1,0.1))
plot(land3, legend = F, axes = F, box = F, col = c("#FFD700", "#548B54"))
plot(land1, legend = F, axes = F, box = F, col = c("#FFD700", "#548B54"))
plot(land4, legend = F, axes = F, box = F, col = c("#FFD700", "#548B54"))
plot(land2, legend = F, axes = F, box = F, col = c("#FFD700", "#548B54"))


######################
### show CTI trend ###
######################

png("../CTI_trend.png", width = 7, height = 4, res = 1200, units = "in")
par(mfrow=c(1,2))

m.FIN <- lmer(cti ~ as.factor(Year) + (1|gridCell50/Site), 
              data = butterflies.data.presence %>% dplyr:::filter(Scale == 50000, country == "FIN"))
m.FIN2 <- lmer(cti ~ Year + (1|gridCell50/Site), 
               data = butterflies.data.presence %>% dplyr:::filter(Scale == 50000, country == "FIN"))

plot(emmean ~ Year, data = emmeans(m.FIN, ~ Year), type = "o", pch = 16, ylab = "Community temperature index")
abline(fixef(m.FIN2)[1], fixef(m.FIN2)[2])

m.NL <- lmer(cti ~ as.factor(Year) + (1|gridCell50/Site), 
             data = butterflies.data.presence %>% dplyr:::filter(Scale == 50000, country == "NL"))
m.NL2 <- lmer(cti ~ Year + (1|gridCell50/Site), 
              data = butterflies.data.presence %>% dplyr:::filter(Scale == 50000, country == "NL"))

plot(emmean ~ Year, data = emmeans(m.NL, ~ Year), type = "o", pch = 16, ylab = "Community temperature index")
abline(fixef(m.NL2)[1], fixef(m.NL2)[2])

dev.off()


###########################
### show classification ###
###########################
library(visreg)

butterflies <- as.tbl(fread("../Data/butterflies_occ.csv"))

par(mfrow=c(2,2))

dat.temp <- butterflies %>% dplyr::filter(Scale == 50000, Site == "1008_NL", Year > 1991)
dat.temp2 <- dat.temp %>% dplyr::filter(Species == "Pararge aegeria")

m <- glm(n ~ Year, family = binomial, dat = dat.temp2)
visreg(m, xvar = "Year", scale = "response", ylim = c(0,1), band = T, rug = F, main = "Colonisation")
points(n ~ Year, dat = dat.temp2)

pred <- predict(m, type = "response", newdata = data.frame(Year = c(min(dat.temp2$Year), max(dat.temp2$Year))))
abline(h = pred[[1]], col = "red")
abline(h = pred[[2]], col = "red")

abline(h = .8, col = "blue", lty = 2)
abline(h = .2, col = "blue", lty = 2)



dat.temp <- butterflies %>% dplyr::filter(Scale == 50000, Site == "101_NL", Year > 1991)
dat.temp2 <- dat.temp %>% dplyr::filter(Species == "Polyommatus icarus")

m <- glm(n ~ Year, family = binomial, dat = dat.temp2)
visreg(m, xvar = "Year", scale = "response", ylim = c(0,1), band = T, rug = F, main = "Extinction")
points(n ~ Year, dat = dat.temp2)

pred <- predict(m, type = "response", newdata = data.frame(Year = c(min(dat.temp2$Year), max(dat.temp2$Year))))
abline(h = pred[[1]], col = "red")
abline(h = pred[[2]], col = "red")

abline(h = .8, col = "blue", lty = 2)
abline(h = .2, col = "blue", lty = 2)



dat.temp <- butterflies %>% dplyr::filter(Scale == 50000, Site == "33_FIN", Year > 1991)
dat.temp2 <- dat.temp %>% dplyr::filter(Species == "Argynnis aglaja")

m <- glm(n ~ Year, family = binomial, dat = dat.temp2)
visreg(m, xvar = "Year", scale = "response", ylim = c(0,1), band = T, rug = F, main = "No colonisation")
points(n ~ Year, dat = dat.temp2)

pred <- predict(m, type = "response", newdata = data.frame(Year = c(min(dat.temp2$Year), max(dat.temp2$Year))))
abline(h = pred[[1]], col = "red")
abline(h = pred[[2]], col = "red")

abline(h = .8, col = "blue", lty = 2)
abline(h = .2, col = "blue", lty = 2)



dat.temp <- butterflies %>% dplyr::filter(Scale == 50000, Site == "1008_NL", Year > 1991)
dat.temp2 <- dat.temp %>% dplyr::filter(Species == "Pieris napi")

m <- glm(n ~ Year, family = binomial, dat = dat.temp2)
visreg(m, xvar = "Year", scale = "response", ylim = c(0,1), band = T, rug = F, main = "Persistence")
points(n ~ Year, dat = dat.temp2)

pred <- predict(m, type = "response", newdata = data.frame(Year = c(min(dat.temp2$Year), max(dat.temp2$Year))))
abline(h = pred[[1]], col = "red")
abline(h = pred[[2]], col = "red")

abline(h = .8, col = "blue", lty = 2)
abline(h = .2, col = "blue", lty = 2)



