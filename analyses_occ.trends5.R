library(data.table)
library(visreg)
library(MuMIn)
library(lmerTest)
library(doFuture)
library(polypoly)
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
plan(multiprocess, workers = 2)

res <- foreach(i = c(1:4)/10, j = rev(c(6:9)/10), .combine = rbind.data.frame) %:% foreach(k = unique(butterflies$Scale), .combine = rbind.data.frame) %dopar%{
  
  dat.Col <- sp_site_Col(dat.occ.trend, butterflies, i, j, 50000)
  dat.Col <- left_join(dat.Col, butterflies %>% dplyr::select(3,10:16,23) %>% group_by(Species) %>% summarise_all(first), by = "Species")
  dat.Col <- left_join(dat.Col, butterflies %>% dplyr::select(2,5,6,8,16:20) %>% group_by(Site, Habitat, Scale) %>% summarise_all(first), by = c("Site", "Habitat"))
  dat.Col <- stdize(dat.Col, prefix = F, omit.cols = c("Extinction", "Colonisation", "Scale", "pred.then", "pred.now", "nYear"))
  
  m.col <- glmer(Colonisation ~ STI_rel * PC3 + 
                   STI_rel * PC4 + 
                   STI_rel * poly_rescale(poly(PC1,2),2) + 
                   poly_rescale(poly(PC1,2),2) * PLAND * CLUMPY + 
                   Habitat + X*Y + 
                   (1|gridCell50/Site) + (1|Species), 
                 family = binomial,
                 data = subset(dat.Col, dat.Col$Scale == k),
                 nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                               optCtrl=list(maxfun=1e10),
                                               calc.derivs = FALSE), verbose = F, na.action = na.fail)
  
  dat.Ext <- sp_site_Ext(dat.occ.trend, i, j)
  dat.Ext <- left_join(dat.Ext, butterflies %>% dplyr::select(2,3,5,6,8,10:20,23) %>% group_by(Species, Site, Scale) %>% summarise_all(first))
  dat.Ext <- stdize(dat.Ext, prefix = F, omit.cols = c("Extinction", "Colonisation", "Scale", "pred.then", "pred.now", "nYear"))
  
  m.ext <- glmer(Extinction ~  STI_rel * PC3 + 
                   STI_rel * PC4 + 
                   STI_rel * poly_rescale(poly(PC1,2),2) + 
                   poly_rescale(poly(PC1,2),2) * PLAND * CLUMPY + 
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

ggplot(res %>% dplyr::filter(grepl("STI_rel|PC1|PC3|PC4|%SNH|Clumpiness", Variables)), 
       aes(x = Variables, y = estimate, color = Scale, alpha = Significancy)) + 
  geom_point(position = position_dodge(.8)) + facet_grid(Process ~ Range) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(.8)) + 
  coord_flip() +  scale_color_brewer(palette="YlOrRd")


### set scale at 20km and range at 0.2 - 0.8 ###
dat.Ext <- sp_site_colExt(dat.occ.trend, .2, .8)
dat.Ext <- left_join(dat.Ext, butterflies %>% dplyr::select(2,3,5,6,8,10:20,23) %>% group_by(Species, Site, Scale) %>% summarise_all(first))
dat.Ext <- stdize(dat.ColExt, prefix = F, omit.cols = c("Extinction", "Colonisation", "Scale", "pred.then", "pred.now", "nYear"))
# scaleList.Ext <- list(scale = attr(dat.ColExt, "scaled:scale"),
#                          center = attr(dat.ColExt, "scaled:center"))


dat.Col <- sp_site_Col(dat.occ.trend, butterflies, .2, .8, 50000)
dat.Col <- left_join(dat.Col, butterflies %>% dplyr::select(3,10:16,23) %>% group_by(Species) %>% summarise_all(first), by = "Species")
dat.Col <- left_join(dat.Col, butterflies %>% dplyr::select(2,5,6,8,16:20) %>% group_by(Site, Habitat, Scale) %>% summarise_all(first), by = c("Site", "Habitat"))
dat.Col <- stdize(dat.Col, prefix = F, omit.cols = c("Extinction", "Colonisation", "Scale", "pred.then", "pred.now", "nYear"))
# scaleList.Col <- list(scale = attr(dat.ColExt, "scaled:scale"),
#                       center = attr(dat.ColExt, "scaled:center"))


##################
## colonisation ##
##################
dat.col.temp <- subset(dat.Col, dat.Col$Scale == 30000)

m.col <- glmer(Colonisation ~ STI_rel * PC3 + 
                 STI_rel * PC4 + 
                 STI_rel * poly_rescale(poly(PC1,2),2) + 
                 poly_rescale(poly(PC1,2),2) * PLAND * CLUMPY + 
                 Habitat + X*Y + 
                 (1|gridCell50/Site) + (1|Species), 
               family = binomial,
               data = dat.col.temp,
               nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                             optCtrl=list(maxfun=1e10),
                                             calc.derivs = FALSE), verbose = F, na.action = na.fail)
summary(m.col)

##################
### Extinction ###
##################
dat.ext.temp <- subset(dat.Ext, dat.Ext$Scale == 30000)

m.ext <- glmer(Extinction ~ STI_rel * PC3 + 
                 STI_rel * PC4 + 
                 STI_rel * PC1 + 
                 PC1 * PLAND * CLUMPY + 
                 Habitat + X*Y + 
                 (1|gridCell50/Site) + (1|Species), 
               family = binomial,
               data = dat.ext.temp,
               nAGQ=0,control = glmerControl(optimizer = "nloptwrap", 
                                             optCtrl=list(maxfun=1e10),
                                             calc.derivs = FALSE), verbose = F)
summary(m.ext)

###########################
### show classification ###
###########################
dat.temp <- butterflies %>% dplyr::filter(Scale == 20000, Site == unique(butterflies$Site)[[600]])

par(mfrow=c(1,2))

dat.temp2 <- dat.temp %>% dplyr::filter(Species == unique(dat.temp$Species)[[15]])

m <- glm(n ~ Year, family = binomial, dat = dat.temp2)
visreg(m, xvar = "Year", scale = "response", ylim = c(0,1), band = T, rug = F)
points(n ~ Year, dat = dat.temp2)

pred <- predict(m, type = "response", newdata = data.frame(Year = c(min(dat.temp2$Year), max(dat.temp2$Year))))
abline(h = pred[[1]], col = "red")
abline(h = pred[[2]], col = "red")

abline(h = .8, col = "blue", lty = 2)
abline(h = .2, col = "blue", lty = 2)


dat.temp2 <- dat.temp %>% dplyr::filter(Species == unique(dat.temp$Species)[[8]])

m <- glm(n ~ Year, family = binomial, dat = dat.temp2)
visreg(m, xvar = "Year", scale = "response", ylim = c(0,1), band = T, rug = F)
points(n ~ Year, dat = dat.temp2)

pred <- predict(m, type = "response", newdata = data.frame(Year = c(min(dat.temp2$Year), max(dat.temp2$Year))))
abline(h = pred[[1]], col = "red")
abline(h = pred[[2]], col = "red")

abline(h = .8, col = "blue", lty = 2)
abline(h = .2, col = "blue", lty = 2)


###########################
## remove unused objects ##
###########################


gdata:::keep(butterflies, dat.occ.trend,
             sure = T)

