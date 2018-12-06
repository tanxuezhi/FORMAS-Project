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

butterflies <- as.tbl(fread("../Data/butterflies_occ.csv"))

dat.occ.trend <- sp_site_occupancy_trend.predict(butterflies %>% group_by(Species, Site, Year) %>% summarise_all(first))
write_csv(dat.occ.trend,"../sp_site_trend_pred.csv")

###############################
### drivers of colonisation ###
###     and extinction at   ###
###   different scales and  ###
###  for various thresholds ###
###############################

dat.occ.trend <- read_csv("../sp_site_trend_pred.csv")



## summary table
table_ColExt <- foreach(i = c(1:4)/10, j = rev(c(6:9)/10), .combine = rbind.data.frame) %do%{
  
  dat <- sp_site_occTrend(left_join(dat.occ.trend, butterflies %>% dplyr::select(2,8) %>% group_by(Site) %>% summarise_all(first)), i, j)
  dat.Col <- dat %>% mutate(Colonisation = ifelse(Trend == "Colonisation", 1, ifelse(Trend == "No colonisation", 0, NA)))
  dat.Ext <- dat %>% mutate(Extinction = ifelse(Trend == "Extinction", 1, ifelse(Trend == "Persistence", 0, NA)))
  cbind.data.frame(thresholds = paste(i, "-", j), 
                   nb.col.YES = table(dat.Col$Colonisation)[2], nb.col.NO = table(dat.Col$Colonisation)[1],
                   nb.ext.YES = table(dat.Ext$Extinction)[2], nb.ext.NO = table(dat.Ext$Extinction)[1])
}


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
                                    confint(m.col, method = "boot")[-c(1:3),]) %>% rownames_to_column("Variables"),
                   cbind.data.frame(Scale = k, 
                                    Process = "Extinction",
                                    Range = paste(i,"-", j), 
                                    estimate = fixef(m.ext), 
                                    confint(m.ext, method = "boot")[-c(1:3),]) %>% rownames_to_column("Variables"))
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
res$Variables <- gsub("PC4", "Resource specialisation", res$Variables)
res$Variables <- gsub("PLAND", "Proportion of SNH", res$Variables)
res$Variables <- gsub("CLUMPY", "Clumpiness of SNH", res$Variables)
res$Variables <- gsub(":", " × ", res$Variables)
res$Variables <- gsub("Proportion of SNH × Clumpiness of SNH", "Proportion × Clumpiness of SNH", res$Variables)

res$Variables <- gsub("STI × Proportion × Clumpiness of SNH", "STI × Proportion ×\nClumpiness of SNH", res$Variables)
res$Variables <- gsub("Spatial use × Proportion × Clumpiness of SNH", "Spatial use × Proportion ×\nClumpiness of SNH", res$Variables)


res <- res %>% mutate(Significancy = ifelse(lwr < 0 & upr < 0 | lwr > 0 & upr > 0, "Significant", "Not significant"))

res$Scale <- factor(res$Scale, levels = rev(paste(c(1,3,5,10,20,30,50), "km")))

res$Variables <- factor(res$Variables, levels = c("(Intercept)", "X", "Y", "X × Y", "HabitatGeneralist", "HabitatOpen",
                                                  "STI", "Spatial use", "Development rate", "Resource specialisation", 
                                                  "STI × Spatial use", "STI × Development rate", "STI × Resource specialisation", 
                                                  "Proportion of SNH", "Clumpiness of SNH", "Proportion × Clumpiness of SNH",
                                                  "Spatial use × Proportion of SNH", "Spatial use × Clumpiness of SNH", 
                                                  "Spatial use × Proportion ×\nClumpiness of SNH",
                                                  "STI × Proportion of SNH", "STI × Clumpiness of SNH", "STI × Proportion ×\nClumpiness of SNH"))



res <- res %>% arrange(Variables, Scale, Process) %>% mutate(Group = rep(c(rep("Various", 6), 
                                                                           rep("Traits", 7),
                                                                           rep("Fragmentation", 3),
                                                                           rep("Fragmentation × traits", 6)), each = 56))

res$Group <- factor(res$Group, levels = c("Various", "Traits", "Fragmentation", "Fragmentation × traits"))

### plot ###

ggplot(res %>% dplyr::filter(grepl("STI|Spatial use|Development rate|Resource specialisation|Proportion|Clumpiness", Variables), 
                             Range == "0.1 - 0.9") %>% 
         mutate(Variables = factor(Variables, levels = rev(levels(Variables)))), 
       aes(x = Variables, y = estimate, color = Scale)) + 
  geom_hline(yintercept = 0, lty = 2, color = "black") +
  geom_point(position = position_dodge(.8)) + facet_grid(Group ~ Process, scales = "free_y", space = "free", switch = "y") +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(.8)) + 
  coord_flip() + 
  scale_y_continuous("Coefficients +/- 95% CI") + 
  scale_x_discrete("") + 
  scale_color_manual(values = colorRampPalette(c("black", "grey"))(7)) + 
  theme_classic() + 
  guides(color = guide_legend(reverse = TRUE)) +
  theme(panel.border=element_blank(),
              strip.text=element_text(size=12, colour="black"),
              strip.background=element_rect(colour="white", 
                                            fill="white"))

cowplot::ggsave("../plot_modelsCoef.svg", width = 10, height = 9)
cowplot::ggsave("../plot_modelsCoef.pdf", width = 10, height = 9)


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
  dplyr::filter(Scale == 20000)
dat.Col$Habitat <- as.factor(dat.Col$Habitat)
dat.Ext <- dat %>% mutate(Extinction = ifelse(Trend == "Extinction", 1, ifelse(Trend == "Persistence", 0, NA))) %>%
  dplyr::filter(Scale == 20000) %>% mutate(Habitat = as.factor(Habitat))
dat.Ext$Habitat <- as.factor(dat.Ext$Habitat)

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
visreg(m.col, xvar = "PC1", by  ="STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["STI_rel"] + scaleList$center["STI_rel"]}, ylim = c(0,.35),
       xlab = "Spatial use", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

visreg(m.col, xvar = "PC4", by  ="STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["STI_rel"] + scaleList$center["STI_rel"]}, ylim = c(0,.35),
       xlab = "Food specialisation", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

plot(0,0,axes = F, ann = F, col = "white")
legend(x = -1, y = 0, title = "STI", legend = c("Low", "High"), lty = 1, col = c("#1C86EE", "#EEC900"), cex = 1.2, lwd = 1.2, bty = "n")

#ext

visreg(m.ext, xvar = "PC1", by = "STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["STI_rel"] + scaleList$center["STI_rel"]},ylim = c(0,.35),
       xlab = "Spatial use", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

visreg(m.ext, xvar = "PC3", by  ="STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["STI_rel"] + scaleList$center["STI_rel"]}, ylim = c(0,.35),
       xlab = "Resource specialisation", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

visreg(m.ext, xvar = "PC4", by  ="STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["PC4"] + scaleList$center["PC4"]}, ylim = c(0,.35),
       xlab = "Development rate", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

dev.off()

# landscape x traits
rm(butterflies); rm(dat); rm(dat.occ.trend); gc()


l1 <- predict_raster3(m.col, xvar = "STI_rel", yvar = "PC1", zvar = "PLAND")
l2 <- predict_raster3(m.ext, xvar = "STI_rel", yvar = "PC1", zvar = "PLAND")
l3 <- predict_raster3(m.col, xvar = "STI_rel", yvar = "PC1", zvar = "CLUMPY")
l4 <- predict_raster3(m.ext, xvar = "STI_rel", yvar = "PC1", zvar = "CLUMPY")


landTraitEffect <- bind_rows(
  "Proportion of SNH" = bind_rows(Colonisation = l1, Extinction = l2, .id = "Process"),
  "Clumpiness of SNH" = bind_rows(Colonisation = l3, Extinction = l4, .id = "Process"),
  .id = "Land_var")

landTraitEffect$Land_var <- as.factor(landTraitEffect$Land_var)

p1 <- ggplot(landTraitEffect %>% filter(Process == "Colonisation", Land_var != "Clumpiness of SNH"), aes(x= STI_rel, y = PC1, fill = estimate)) +
  geom_raster() + facet_grid(Land_var~., scales = "free", drop = F) +
  scale_fill_gradientn(name = "Effect on \ncolonisation\nprobability\n",
                       colours=c("white","blue")) +
  scale_x_continuous("Species temperature index") + scale_y_continuous("Spatial use") +
  theme_classic() +
  theme(panel.border=element_blank(),
        strip.text=element_text(size=12, colour="black"),
        strip.background=element_rect(colour="white", 
                                      fill="white"))

p2 <- ggplot(landTraitEffect %>% filter(Process == "Extinction"), aes(x= STI_rel, y = PC1, fill = estimate)) +
  geom_raster() + facet_grid(Land_var~., scales = "free") +
  scale_fill_gradientn(name = "Effect on \nextinction\nprobability\n",
                       colours=c("yellow","red")) +
  scale_x_continuous("Species temperature index") + scale_y_continuous("Spatial use") +
  theme_classic() +
  theme(panel.border=element_blank(),
        strip.text=element_text(size=12, colour="black"),
        strip.background=element_rect(colour="white", 
                                      fill="white"))

p3 <- cowplot::plot_grid(p1, p2, labels = c("Colonisation", "Extinction"))

cowplot::ggsave("../plot_traitsxLand3.pdf", p3, width = 8, height = 4.5, dpi = 300)


# alt
landTraitEffect <- bind_rows(
  "Low" = bind_rows(Colonisation = predict_raster2(m.col, xvar = "PLAND", yvar = "CLUMPY", cond = list(PC1 = -1)),
                              Extinction = predict_raster2(m.ext, xvar = "PLAND", yvar = "CLUMPY", cond = list(PC1 = -1)), .id = "Process"),
  "High" = bind_rows(Colonisation = predict_raster2(m.col, xvar = "PLAND", yvar = "CLUMPY", cond = list(PC1 = 1)),
                               Extinction = predict_raster2(m.ext, xvar = "PLAND", yvar = "CLUMPY", cond = list(PC1 = 1)), .id = "Process"),
  .id = "Value")

landTraitEffect2 <- bind_rows(
  "Low" = bind_rows(Colonisation = predict_raster2(m.col, xvar = "PLAND", yvar = "CLUMPY", cond = list(STI_rel = -1)),
                              Extinction = predict_raster2(m.ext, xvar = "PLAND", yvar = "CLUMPY", cond = list(STI_rel = -1)), .id = "Process"),
  "High" = bind_rows(Colonisation = predict_raster2(m.col, xvar = "PLAND", yvar = "CLUMPY", cond = list(STI_rel = 1)),
                               Extinction = predict_raster2(m.ext, xvar = "PLAND", yvar = "CLUMPY", cond = list(STI_rel = 1)), .id = "Process"),
  .id = "Value")

landTraitEffectfull <- bind_rows("Spatial use" = landTraitEffect, "STI" = landTraitEffect2, .id = "Trait")
landTraitEffectfull$Value <- factor(landTraitEffectfull$Value, levels = c("Low", "High"))


ggplot(landTraitEffectfull %>% filter(Process == "Extinction"), aes(x= PLAND, y = CLUMPY, fill = log10(pred))) +
  geom_raster() + facet_grid(Value~Trait) +
  scale_fill_gradientn(name = "Extinction\nprobability\n",
                       colours=c("yellow","red"),
                       limits = c(-5.25, 0),
                       labels = function(x){10^x}) +
  scale_x_continuous("Proportion of SNH") + scale_y_continuous("Clumpiness of SNH") +
  theme_classic() +
  theme(panel.border=element_blank(),
                       strip.text=element_text(size=12, colour="black"),
                       strip.background=element_rect(colour="white", 
                                                     fill="white"))

ggplot(landTraitEffectfull %>% filter(Process == "Colonisation"), aes(x= PLAND, y = CLUMPY, fill = log10(pred))) +
  geom_raster() + facet_grid(Value~Trait) +
  scale_fill_gradientn(name = "Colonisation\nprobability\n",
                       colours=c("white","dark blue"),
                       limits = c(-2.6, 0),
                       labels = function(x){10^x}) +
  scale_x_continuous("Proportion of SNH") + scale_y_continuous("Clumpiness of SNH") +
  theme_classic() + geom_contour(aes(z = pred)) + 
  theme(panel.border=element_blank(),
        strip.text=element_text(size=12, colour="black"),
        strip.background=element_rect(colour="white", 
                                      fill="white"))




svg("../plot_traitsxLand2.svg", width = 7, height = 2.5)


par(mfrow=c(2,4), mar = c(4, 4, 2.5, 2))
visreg(m.col, xvar = "PLAND", by  ="PC1", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["PLAND"] + scaleList$center["PLAND"]}, ylim = c(0,.35),
       xlab = "Proportion of SNH", ylab = "Colonisation probability", legend = T,
       line.par=list(col=c("#8DB6CD", "#7A378B")))

visreg(m.col, xvar = "CLUMPY", by  ="PC1", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"]}, ylim = c(0,.35),
       xlab = "Clumpiness of SNH", ylab = "Colonisation probability", legend = T, 
       line.par = list(col = c("#8DB6CD", "#7A378B")))

visreg(m.ext, xvar = "PLAND", by  ="PC1", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["PLAND"] + scaleList$center["PLAND"]}, ylim = c(0,.35),
       xlab = "Proportion of SNH", ylab = "Extinction probability", legend = T,
       line.par=list(col=c("#8DB6CD", "#7A378B")))

visreg(m.ext, xvar = "CLUMPY", by  ="PC1", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"]}, ylim = c(0,.35),
       xlab = "Clumpiness of SNH", ylab = "Extinction probability", legend = T, 
       line.par = list(col = c("#8DB6CD", "#7A378B")))



visreg(m.col, xvar = "PLAND", by  ="STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["PLAND"] + scaleList$center["PLAND"]}, ylim = c(0,.35),
       xlab = "Proportion of SNH", ylab = "Colonisation probability", legend = T,
       line.par=list(col=c("#8DB6CD", "#7A378B")))

visreg(m.col, xvar = "CLUMPY", by  ="STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"]}, ylim = c(0,.35),
       xlab = "Clumpiness of SNH", ylab = "Colonisation probability", legend = T, 
       line.par = list(col = c("#8DB6CD", "#7A378B")))

visreg(m.ext, xvar = "PLAND", by  ="STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["PLAND"] + scaleList$center["PLAND"]}, ylim = c(0,.35),
       xlab = "Proportion of SNH", ylab = "Extinction probability", legend = T,
       line.par=list(col=c("#8DB6CD", "#7A378B")))

visreg(m.ext, xvar = "CLUMPY", by  ="STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = F,
       xtrans = function(x){x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"]}, ylim = c(0,.35),
       xlab = "Clumpiness of SNH", ylab = "Extinction probability", legend = T, 
       line.par = list(col = c("#8DB6CD", "#7A378B")))


dev.off()

###########################
## remove unused objects ##
###########################


gdata:::keep(butterflies, dat.occ.trend,
             sure = T)

