library(data.table)
library(visreg)
library(MuMIn)
library(lmerTest)
library(doFuture)
library(patchwork)
library(tidyverse)
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
  
  dat <- sp_site_occTrend(left_join(dat.occ.trend, butterflies %>% dplyr::select(2,8) %>% group_by(Site) %>% summarise_all(first)), i, j)
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
res$Variables <- gsub("PC1", "Spatial use", res$Variables)
res$Variables <- gsub("PC3", "Development rate", res$Variables)
res$Variables <- gsub("PC4", "Food specialisation", res$Variables)
res$Variables <- gsub("PLAND", "%SNH", res$Variables)
res$Variables <- gsub("CLUMPY", "Clumpiness", res$Variables)
res$Variables <- gsub(":", " x ", res$Variables)
res <- res %>% mutate(Significancy = ifelse(lwr < 0 & upr < 0 | lwr > 0 & upr > 0, "Significant", "Not significant"))

res$Scale <- factor(res$Scale, levels = paste(c(1,3,5,10,20,30,50), "km"))

res$Variables <- factor(res$Variables, levels = c("(Intercept)", "X", "Y", "X x Y", "HabitatGeneralist", "HabitatOpen",
                                                  "STI", "Spatial use", "Development rate", "Food specialisation", 
                                                  "STI x Spatial use", "STI x Development rate", "STI x Food specialisation", 
                                                  "%SNH", "Clumpiness", "%SNH x Clumpiness",
                                                  "Spatial use x %SNH", "Spatial use x Clumpiness", 
                                                  "Spatial use x %SNH x Clumpiness",
                                                  "STI x %SNH", "STI x Clumpiness", "STI x %SNH x Clumpiness"))
res <- res %>% arrange(Variables, Scale, Process) %>% mutate(Group = rep(c(rep("Various", 6), 
                                                                           rep("Traits", 7),
                                                                           rep("Landscape", 3),
                                                                           rep("Landscape x traits", 6)), each = 56))

res$Group <- factor(res$Group, levels = c("Various", "Traits", "Landscape", "Landscape x traits"))

### plot ###

ggplot(res %>% dplyr::filter(grepl("STI|Spatial use|Development rate|Food specialisation|%SNH|Clumpiness", Variables), 
                             Range == "0.2 - 0.8") %>% 
         mutate(Variables = factor(Variables, levels = rev(levels(Variables)))), 
       aes(x = Variables, y = estimate, color = Scale, alpha = Significancy)) + 
  geom_point(position = position_dodge(.8)) + facet_grid(Group ~ Process, scales = "free_y", space = "free", switch = "y") +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(.8)) + 
  geom_hline(yintercept = 0, lty = 2, color = "grey") +
  coord_flip() + 
  scale_color_manual(values = c('#eff3ff','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#084594')) + 
  scale_alpha_discrete(guide=FALSE, range = c(.2, 1)) + 
  theme_bw()

ggsave("../plot_modelsCoef.svg", width = 8, height = 6)

###############################
### drivers of colonisation ###
###     and extinction at   ###
###      20km scale and     ###
###    .2 - .8 threshold    ###
###############################


### extract data at scale = 20km and range = 0.2 - 0.8 ###
dat <- sp_site_occTrend(left_join(dat.occ.trend, butterflies %>% dplyr::select(2,8) %>% group_by(Site) %>% summarise_all(first)), .2, .8)
table(dat$Trend)

dat <- left_join(dat, butterflies %>% dplyr::select(2,3,5,6,8,10:20,23) %>% group_by(Species, Site, Scale) %>% summarise_all(first))
dat <- stdize(dat, prefix = F, omit.cols = c("Trend", "Scale", "pred.then", "pred.now", "nYear"))

dat.Col <- dat %>% mutate(Colonisation = ifelse(Trend == "Colonisation", 1, ifelse(Trend == "No colonisation", 0, NA))) %>%
  dplyr::filter(Scale == 30000)
dat.Ext <- dat %>% mutate(Extinction = ifelse(Trend == "Extinction", 1, ifelse(Trend == "Persistence", 0, NA))) %>%
  dplyr::filter(Scale == 30000)

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
# landEffect <- bind_rows(Colonisation = predict_raster2(m.col, xvar = "PLAND", yvar = "CLUMPY"),
#                         Extinction = predict_raster2(m.ext, xvar = "PLAND", yvar = "CLUMPY"), .id = "Process")
# 
# ggplot(landEffect, aes(x= PLAND, y = CLUMPY, fill = pred)) + geom_raster() + facet_wrap(~Process) +
#   scale_fill_gradientn(name = "Colonisation /\nextinction\nprobability\n",
#                        breaks =  c(min(landEffect$pred),max(landEffect$pred)),
#                        labels = c(round(min(landEffect$pred), 3),round(max(landEffect$pred), 3)),
#                        limits = c(min(landEffect$pred),max(landEffect$pred)),
#                        colours=c("yellow","red")) +
#   scale_x_continuous("% Semi-natural habitat") + scale_y_continuous("Habitat clumpiness")


# png("../plot_landscape.png", width = 3.5, height = 3.5, units = "in", res = 300)
# 
# par(mfrow=c(1,1), mar = c(4, 4, 1, 2))
# visreg(m.ext, xvar = "CLUMPY", by = "PLAND", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = T,
#        xtrans = function(x){x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"]}, ylim = c(0,.3),
#        xlab = "Clumpiness of SNH", legend = F, fill.par = list(col = c("#8DB6CD60", "#7A378B60")), line.par = list(col = c("#8DB6CD", "#7A378B")))
# legend("topright", title = "Proportion of SNH", legend = c("Low", "High"), lty = 1, col = c("#8DB6CD", "#7A378B"), lwd = 1.2, bty = "n")
# dev.off()



# traits only
svg("../plot_traits2.svg", width = 6, height = 4)

par(mfrow=c(2,3), mar = c(4, 4, 2.5, 2))
#col
visreg(m.col, xvar = "PC1", by  ="STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = T,
       xtrans = function(x){x * scaleList$scale["STI_rel"] + scaleList$center["STI_rel"]}, ylim = c(0,.4),
       xlab = "Spatial use", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

visreg(m.col, xvar = "PC4", by  ="STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = T,
       xtrans = function(x){x * scaleList$scale["STI_rel"] + scaleList$center["STI_rel"]}, ylim = c(0,.4),
       xlab = "Food specialisation", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

plot(0,0,axes = F, ann = F, col = "white")
legend(x = -1, y = 0, title = "STI", legend = c("Low", "High"), lty = 1, col = c("#1C86EE", "#EEC900"), cex = 1.2, lwd = 1.2, bty = "n")

#ext

visreg(m.ext, xvar = "PC1", by = "STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = T,
       xtrans = function(x){x * scaleList$scale["STI_rel"] + scaleList$center["STI_rel"]},ylim = c(0,.15),
       xlab = "Spatial use", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

visreg(m.ext, xvar = "PC3", by  ="STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = T,
       xtrans = function(x){x * scaleList$scale["STI_rel"] + scaleList$center["STI_rel"]}, ylim = c(0,.15),
       xlab = "Food specialisation", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

visreg(m.ext, xvar = "PC4", by  ="STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = T,
       xtrans = function(x){x * scaleList$scale["PC4"] + scaleList$center["PC4"]}, ylim = c(0,.15),
       xlab = "Development rate", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

dev.off()

# landscape x traits

# alt
# landTraitEffect <- bind_rows(
#   "Low dispersal" = bind_rows(Colonisation = predict_raster2(m.col, xvar = "PLAND", yvar = "CLUMPY", cond = list(PC1 = -1)),
#                               Extinction = predict_raster2(m.ext, xvar = "PLAND", yvar = "CLUMPY", cond = list(PC1 = -1)), .id = "Process"),
#   "High dispersal" = bind_rows(Colonisation = predict_raster2(m.col, xvar = "PLAND", yvar = "CLUMPY", cond = list(PC1 = 1)),
#                                Extinction = predict_raster2(m.ext, xvar = "PLAND", yvar = "CLUMPY", cond = list(PC1 = 1)), .id = "Process"),
#   .id = "PC1")
# 
# landTraitEffect2 <- bind_rows(
#   "Low STI" = bind_rows(Colonisation = predict_raster2(m.col, xvar = "PLAND", yvar = "CLUMPY", cond = list(STI_rel = -1)),
#                               Extinction = predict_raster2(m.ext, xvar = "PLAND", yvar = "CLUMPY", cond = list(STI_rel = -1)), .id = "Process"),
#   "High STI" = bind_rows(Colonisation = predict_raster2(m.col, xvar = "PLAND", yvar = "CLUMPY", cond = list(STI_rel = 1)),
#                                Extinction = predict_raster2(m.ext, xvar = "PLAND", yvar = "CLUMPY", cond = list(STI_rel = 1)), .id = "Process"),
#   .id = "STI_rel")
# 
# 
# ggplot(landTraitEffect %>% filter(Process == "Extinction"), aes(x= PLAND, y = CLUMPY, fill = pred)) + 
#   geom_raster() + facet_grid(.~PC1) +
#   scale_fill_gradientn(name = "Extinction\nprobability\n",
#                        breaks =  c(min(landTraitEffect$pred),max(landTraitEffect$pred)),
#                        labels = c(round(min(landTraitEffect$pred), 3),round(max(landTraitEffect$pred), 3)),
#                        limits = c(min(landTraitEffect$pred),max(landTraitEffect$pred)),
#                        colours=c("yellow","red")) +
#   scale_x_continuous("% Semi-natural habitat") + scale_y_continuous("Habitat clumpiness") +
#   
# ggplot(landTraitEffect2 %>% filter(Process == "Extinction"), aes(x= PLAND, y = CLUMPY, fill = pred)) + 
#   geom_raster() + facet_grid(.~STI_rel) +
#   scale_fill_gradientn(name = "Extinction\nprobability\n",
#                        breaks =  c(min(landTraitEffect2$pred),max(landTraitEffect2$pred)),
#                        labels = c(round(min(landTraitEffect2$pred), 3),round(max(landTraitEffect2$pred), 3)),
#                        limits = c(min(landTraitEffect2$pred),max(landTraitEffect2$pred)),
#                        colours=c("yellow","red")) +
#   scale_x_continuous("% Semi-natural habitat") + scale_y_continuous("Habitat clumpiness")


svg("../plot_traitsxLand2.svg", width = 7, height = 2.5)

par(mfrow=c(1,3), mar = c(4, 4, 2.5, 2))
visreg(m.col, xvar = "PLAND", by  ="PC1", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = T,
       xtrans = function(x){x * scaleList$scale["PLAND"] + scaleList$center["PLAND"]}, ylim = c(0,.4),
       xlab = "Proportion of SNH", ylab = "Colonisation probability", legend = T, fill.par=list(col=c("#8DB6CD70", "#7A378B40")),
       line.par=list(col=c("#8DB6CD", "#7A378B")))

visreg(m.ext, xvar = "CLUMPY", by  ="PLAND", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = T,
       xtrans = function(x){x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"]}, ylim = c(0,.5),
       xlab = "Aggregation of SNH", ylab = "Extinction probability", legend = F, fill.par = list(col = c("#7A378B40", "#7A378B40")), line.par = list(col = c("#7A378B", "#7A378B"), lty = c(3,1)), cond = list(PC1 = 1))
par(new = T)
visreg(m.ext, xvar = "CLUMPY", by  ="PLAND", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = T,
       xtrans = function(x){x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"]}, ylim = c(0,.5),
       xlab = "", ylab = "", legend = F, fill.par = list(col = c("#8DB6CD70", "#8DB6CD70")),ann = F, axes = F,
       line.par = list(col = c("#8DB6CD", "#8DB6CD"), lty = c(3,1)), cond = list(PC1 = -1))


visreg(m.ext, xvar = "CLUMPY", by  ="PLAND", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = T,
       xtrans = function(x){x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"]}, ylim = c(0,.5),
       xlab = "Aggregation of SNH", ylab = "Extinction probability", legend = F, fill.par=list(col=c("#1C86EE60", "#1C86EE60")),
       line.par=list(col=c("#1C86EE", "#1C86EE"), lty = c(3,1)), cond = list(STI_rel = -1))
par(new = T)
visreg(m.ext, xvar = "CLUMPY", by  ="PLAND", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = T,
       xtrans = function(x){x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"]}, ylim = c(0,.5),
       xlab = "", ylab = "", legend = F, fill.par=list(col=c("#EEC90060", "#EEC90060")),ann = F, axes = F,
       line.par=list(col=c("#EEC900", "#EEC900"), lty = c(3,1)), cond = list(STI_rel = 1))
  
dev.off()

###########################
## remove unused objects ##
###########################


gdata:::keep(butterflies, dat.occ.trend,
             sure = T)

