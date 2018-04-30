library(data.table)
library(visreg)
library(MuMIn)
library(lmerTest)
library(doFuture)
source("functions.R")

###############################
### extract occupancy trend ###
###   by species and site   ###
###############################

butterflies <- as.tbl(fread("C:/Local Folder (c)/butterflies_occ.csv"))

dat.occ.trend <- sp_site_occupancy_trend.predict(butterflies %>% group_by(Species, Site, Year) %>% summarise_all(first))
write_csv(dat.occ.trend,"../sp_site_trend_pred.csv")

###############################
### drivers of colonisation ###
###     and extinction at   ###
###   different scales and  ###
###  for various thresholds ###
###############################

dat.occ.trend <- read_csv("../sp_site_trend_pred.csv")

registerDoFuture()
options(future.globals.maxSize = Inf)
plan(multiprocess, workers = 3)

res <- foreach(i = c(1:4)/10, j = rev(c(6:9)/10), .combine = rbind.data.frame) %:% foreach(k = unique(butterflies$Scale), .combine = rbind.data.frame) %dopar%{
  
  dat <- sp_site_occTrend(dat.occ.trend, i, j)
  dat <- left_join(dat, butterflies %>% dplyr::select(2,3,5,6,8,10:20,23) %>% group_by(Species, Site, Scale) %>% summarise_all(first))
  dat <- stdize(dat, prefix = F, omit.cols = c("Trend", "Scale", "pred.then", "pred.now", "nYear"))
  
  dat.Col <- dat %>% mutate(Colonisation = ifelse(Trend == "Colonisation", 1, ifelse(Trend == "No colonisation", 0, NA)))
  dat.Ext <- dat %>% mutate(Extinction = ifelse(Trend == "Extinction", 1, ifelse(Trend == "Persistence", 0, NA)))
  
  m.col <- glmer(Colonisation ~ STI_rel * PC3 + 
                   STI_rel * PC4 + 
                   STI_rel * PC1 + 
                   PC1 * PLAND * CLUMPY + 
                   STI_rel * PLAND * CLUMPY + 
                   Habitat + X*Y + 
                   (1|gridCell50/Site) + (1|Species), 
                 family = binomial,
                 data = subset(dat.Col, dat.Col$Scale == k),
                 nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                               optCtrl=list(maxfun=1e10),
                                               calc.derivs = FALSE), verbose = F)
  
  m.ext <- glmer(Extinction ~  STI_rel * PC3 + 
                   STI_rel * PC4 + 
                   STI_rel * PC1 + 
                   PC1 * PLAND * CLUMPY + 
                   STI_rel * PLAND * CLUMPY + 
                   Habitat + X*Y + 
                   (1|gridCell50/Site) + (1|Species), 
                 family = binomial,
                 data = subset(dat.Ext, dat.Ext$Scale == k),
                 nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                               optCtrl=list(maxfun=1e10),
                                               calc.derivs = FALSE), verbose = F)
  
  rbind.data.frame(cbind.data.frame(Scale = k, 
                                    Process = "Colonisation",
                                    Range = paste(i,"-", j), 
                                    estimate = fixef(m.col), 
                                    confint(m.col, method = "Wald")[-c(1:3),]) %>% rownames_to_column("Variables"),
                   cbind.data.frame(Scale = k, 
                                    Process = "Extinction",
                                    Range = paste(i,"-", j), 
                                    estimate = fixef(m.ext), 
                                    confint(m.ext, method = "Wald")[-c(1:3),]) %>% rownames_to_column("Variables"))
}

names(res)[6:7] <- c("lwr", "upr")

write_csv(res, "../coefs_mods_ColExt.csv")

res <- read_csv("../coefs_mods_ColExt.csv")
res$Scale <- as.factor(res$Scale)
levels(res$Scale) <- paste(c(1,3,5,10,20,30,50), "km")

res$Variables <- gsub("poly_rescale\\(poly\\(PC1, 2\\), 2\\)2", "PC1^2", res$Variables)
res$Variables <- gsub("poly_rescale\\(poly\\(PC1, 2\\), 2\\)1", "PC1", res$Variables)
res$Variables <- gsub("STI_rel", "STI", res$Variables)
res$Variables <- gsub("PLAND", "%SNH", res$Variables)
res$Variables <- gsub("CLUMPY", "Clumpiness", res$Variables)
res$Variables <- gsub(":", " x ", res$Variables)
res <- res %>% mutate(Significancy = ifelse(lwr < 0 & upr < 0 | lwr > 0 & upr > 0, "Significant", "Not significant"))

res$Scale <- factor(res$Scale, levels = paste(c(1,3,5,10,20,30,50), "km"))

res$Variables <- factor(res$Variables, levels = c("(Intercept)", "X", "Y", "X x Y", "HabitatGeneralist", "HabitatOpen",
                                                  "STI", "PC1", "PC3", "PC4", 
                                                  "STI x PC1", "STI x PC3", "STI x PC4", 
                                                  "%SNH", "Clumpiness", "%SNH x Clumpiness",
                                                  "PC1 x %SNH", "PC1 x Clumpiness", "PC1 x %SNH x Clumpiness",
                                                  "STI x %SNH", "STI x Clumpiness", "STI x %SNH x Clumpiness"))
res <- res %>% arrange(Variables, Scale, Process) %>% mutate(Group = rep(c(rep("Various", 6), 
                                                                           rep("Traits", 7),
                                                                           rep("Landscape", 3),
                                                                           rep("Landscape x traits", 6)), each = 56))

res$Group <- factor(res$Group, levels = c("Various", "Traits", "Landscape", "Landscape x traits"))

### plot ###

ggplot(res %>% dplyr::filter(grepl("STI|PC1|PC3|PC4|%SNH|Clumpiness", Variables), Range == "0.2 - 0.8") %>% 
         mutate(Variables = factor(Variables, levels = rev(levels(Variables)))), 
       aes(x = Variables, y = estimate, color = Scale, alpha = Significancy)) + 
  geom_point(position = position_dodge(.8)) + facet_grid(Group ~ Process, scales = "free_y", space = "free", switch = "y") +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(.8)) + 
  geom_hline(yintercept = 0, lty = 2, color = "grey") +
  coord_flip() + 
  scale_color_brewer(palette="YlOrRd") + 
  scale_alpha_discrete(guide=FALSE, range = c(.2, 1)) + 
  theme_bw()


###############################
### drivers of colonisation ###
###     and extinction at   ###
###      20km scale and     ###
###    .2 - .8 threshold    ###
###############################


### extract data at scale = 20km and range = 0.2 - 0.8 ###
dat <- sp_site_occTrend(dat.occ.trend, .2, .8)
table(dat$Trend)

dat <- left_join(dat, butterflies %>% dplyr::select(2,3,5,6,8,10:20,23) %>% group_by(Species, Site, Scale) %>% summarise_all(first))
dat <- stdize(dat, prefix = F, omit.cols = c("Trend", "Scale", "pred.then", "pred.now", "nYear"))

dat.Col <- dat %>% mutate(Colonisation = ifelse(Trend == "Colonisation", 1, ifelse(Trend == "No colonisation", 0, NA))) %>%
  dplyr::filter(Scale == 20000)
dat.Ext <- dat %>% mutate(Extinction = ifelse(Trend == "Extinction", 1, ifelse(Trend == "Persistence", 0, NA))) %>%
  dplyr::filter(Scale == 20000)

scaleList <- list(scale = attr(dat, "scaled:scale"),
                  center = attr(dat, "scaled:center"))


##################
## colonisation ##
##################
m.col <- glmer(Colonisation ~ STI_rel * PC3 + 
                 STI_rel * PC4 + 
                 STI_rel * PC1 + 
                 PC1 * PLAND * CLUMPY + 
                 STI_rel * PLAND * CLUMPY + 
                 Habitat + X*Y + 
                 (1|gridCell50/Site) + (1|Species), 
               family = binomial,
               data = dat.Col,
               nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                             optCtrl=list(maxfun=1e10),
                                             calc.derivs = FALSE), verbose = F)
summary(m.col)

##################
### Extinction ###
##################
m.ext <- glmer(Extinction ~ STI_rel * PC3 + 
                 STI_rel * PC4 + 
                 STI_rel * PC1 + 
                 PC1 * PLAND * CLUMPY + 
                 STI_rel * PLAND * CLUMPY + 
                 Habitat + X*Y + 
                 (1|gridCell50/Site) + (1|Species), 
               family = binomial,
               data = dat.Ext,
               nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                             optCtrl=list(maxfun=1e10),
                                             calc.derivs = FALSE), verbose = F)
summary(m.ext)


## plots ##

# landscape only
landEffect <- bind_rows(Colonisation = predict_raster2(m.col, xvar = "PLAND", yvar = "CLUMPY"), 
                        Extinction = predict_raster2(m.ext, xvar = "PLAND", yvar = "CLUMPY"), .id = "Process")

ggplot(landEffect, aes(x= PLAND, y = CLUMPY, fill = pred)) + geom_raster() + facet_wrap(~Process) +
  scale_fill_gradientn(name = "Colonisation /\nextinction\nprobability\n", 
                       breaks =  c(min(landEffect$pred),max(landEffect$pred)),
                       labels = c(round(min(landEffect$pred), 3),round(max(landEffect$pred), 3)),
                       limits = c(min(landEffect$pred),max(landEffect$pred)),
                       colours=c("yellow","red")) +
  scale_x_continuous("% Semi-natural habitat") + scale_y_continuous("Habitat clumpiness")


# traits only
png("../plot.png", width = 4, height = 7, units = "in", res = 300)

par(mfrow=c(3,2))
#col
visreg(m.col, xvar = "STI_rel", by  ="PC4", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["PC4"] + scaleList$center["PC4"]}, ylim = c(0,1),
       xlab = "STI", legend = F)
legend("topleft", title = "PC4", legend = c("Low", "High"), lty = 1, col = c("red", "blue"), bty = "n")

visreg(m.col, xvar = "STI_rel", by  ="PC1", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["PC4"] + scaleList$center["PC4"]}, ylim = c(0,1),
       xlab = "STI", legend = F)
legend("topleft", title = "PC1", legend = c("Low", "High"), lty = 1, col = c("red", "blue"), bty = "n")

#ext
# visreg(m.ext, xvar = "STI_rel", by  ="STI_sd", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
#        xtrans = function(x){x * scaleList$scale["STI_rel"] + scaleList$center["STI_rel"]},ylim = c(0,1),
#        xlab = "STI", legend = F)
# legend("topright", title = "Std. dev STI", legend = c("Low", "High"), lty = 1, col = c("red", "blue"), bty = "n")

visreg(m.ext, xvar = "STI_rel", by  ="PC1", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["STI_rel"] + scaleList$center["STI_rel"]},ylim = c(0,1),
       xlab = "STI", legend = F)
legend("topright", title = "PC1", legend = c("Low", "High"), lty = 1, col = c("red", "blue"), bty = "n")

visreg(m.ext, xvar = "STI_rel", by  ="PC3", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["STI_rel"] + scaleList$center["STI_rel"]}, ylim = c(0,1),
       xlab = "STI", legend = F)
legend("topright", title = "PC3", legend = c("Low", "High"), lty = 1, col = c("red", "blue"), bty = "n")

visreg(m.ext, xvar = "STI_rel", by  ="PC4", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["PC4"] + scaleList$center["PC4"]}, ylim = c(0,1),
       xlab = "STI", legend = F)
legend("topright", title = "PC4", legend = c("Low", "High"), lty = 1, col = c("red", "blue"), bty = "n")

dev.off()

# landscape x traits
png("../plot.png", width = 4, height = 7, units = "in", res = 300)

par(mfrow=c(3,2))
visreg(m.col, xvar = "PLAND", by  ="PC1", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["PLAND"] + scaleList$center["PLAND"]}, ylim = c(0,0.3),
       strip.names = c("Low PC1", "High PC1"))

plot.new()

visreg(m.ext, xvar = "PLAND", by  ="PC1", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["PLAND"] + scaleList$center["PLAND"]}, ylim = c(0,0.3),
       strip.names = c("Low PC1", "High PC1"))

visreg(m.ext, xvar = "CLUMPY", by  ="PC1", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"]}, ylim = c(0,0.3),
       strip.names = c("Low PC1", "High PC1"))

visreg(m.ext, xvar = "CLUMPY", by  ="STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"]}, ylim = c(0,0.3),
       strip.names = c("Low STI", "High STI"))

visreg(m.ext, xvar = "PLAND", by  ="STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["PLAND"] + scaleList$center["PLAND"]}, ylim = c(0,0.3),
       strip.names = c("Low STI", "High STI"))

dev.off()

# alt
landTraitEffect <- bind_rows(
  "Low dispersal" = bind_rows(Colonisation = predict_raster2(m.col, xvar = "PLAND", yvar = "CLUMPY", cond = list(PC1 = -1)), 
                              Extinction = predict_raster2(m.ext, xvar = "PLAND", yvar = "CLUMPY", cond = list(PC1 = -1)), .id = "Process"),
  "High dispersal" = bind_rows(Colonisation = predict_raster2(m.col, xvar = "PLAND", yvar = "CLUMPY", cond = list(PC1 = 1)), 
                               Extinction = predict_raster2(m.ext, xvar = "PLAND", yvar = "CLUMPY", cond = list(PC1 = 1)), .id = "Process"),
  .id = "PC1")

ggplot(landTraitEffect, aes(x= PLAND, y = CLUMPY, fill = pred)) + geom_raster() + facet_grid(PC1~Process) +
  scale_fill_gradientn(name = "Colonisation /\nextinction\nprobability\n", 
                       breaks =  c(min(landTraitEffect$pred),max(landTraitEffect$pred)),
                       labels = c(round(min(landTraitEffect$pred), 3),round(max(landTraitEffect$pred), 3)),
                       limits = c(min(landTraitEffect$pred),max(landTraitEffect$pred)),
                       colours=c("yellow","red")) +
  scale_x_continuous("% Semi-natural habitat") + scale_y_continuous("Habitat clumpiness")


###########################
### show classification ###
###########################
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


# dat.temp <- butterflies %>% dplyr::filter(Scale == 50000, Site == "81_FIN", Year > 1991)
# dat.temp2 <- dat.temp %>% dplyr::filter(Species == "Argynnis aglaja")
# 
# m <- glm(n ~ Year, family = binomial, dat = dat.temp2)
# visreg(m, xvar = "Year", scale = "response", ylim = c(0,1), band = T, rug = F, main = "Undetermined")
# points(n ~ Year, dat = dat.temp2)
# 
# pred <- predict(m, type = "response", newdata = data.frame(Year = c(min(dat.temp2$Year), max(dat.temp2$Year))))
# abline(h = pred[[1]], col = "red")
# abline(h = pred[[2]], col = "red")
# 
# abline(h = .8, col = "blue", lty = 2)
# abline(h = .2, col = "blue", lty = 2)


###########################
## remove unused objects ##
###########################


gdata:::keep(butterflies, dat.occ.trend,
             sure = T)

