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


## summary table
table_ColExt <- foreach(i = c(1:4)/10, j = rev(c(6:9)/10), .combine = rbind.data.frame) %do%{
  
  dat <- sp_site_occTrend(left_join(dat.occ.trend, butterflies %>% dplyr::select(2,8) %>% group_by(Site) %>% summarise_all(first)), i, j)
  dat.Col <- dat %>% mutate(Colonisation = ifelse(Trend == "Colonisation", 1, ifelse(Trend == "No colonisation", 0, NA)))
  dat.Ext <- dat %>% mutate(Extinction = ifelse(Trend == "Extinction", 1, ifelse(Trend == "Persistence", 0, NA)))
  
  cbind.data.frame(thresholds = paste(i, "-", j), rbind.data.frame(
    cbind.data.frame(Type = "Colonisation",
                     dat.Col %>% group_by(Country = ifelse(grepl("_NL", Site), "NL", "FIN"), Colonisation) %>% 
                       summarise(n = n()) %>% rename(yes_no = Colonisation)),
    cbind.data.frame(Type = "Extinction",
                     dat.Ext %>% group_by(Country = ifelse(grepl("_NL", Site), "NL", "FIN"), Extinction) %>% 
                       summarise(n = n())%>% rename(yes_no = Extinction))
  )
  )
}

write.table(table_ColExt, "clipboard")

# summary table per species

table_ColExt_sp <- foreach(i = c(1:4)/10, j = rev(c(6:9)/10), .combine = rbind.data.frame) %do%{
  
  dat <- sp_site_occTrend(left_join(dat.occ.trend, butterflies %>% dplyr::select(2,8) %>% 
                                      group_by(Site) %>% 
                                      summarise_all(first)), i, j)
  dat.Col <- dat %>% mutate(Colonisation = ifelse(Trend == "Colonisation", 1, ifelse(Trend == "No colonisation", 0, NA)),
                            Country = ifelse(grepl("_NL", Site), "NL", "FIN"))
  dat.Ext <- dat %>% mutate(Extinction = ifelse(Trend == "Extinction", 1, ifelse(Trend == "Persistence", 0, NA)),
                            Country = ifelse(grepl("_NL", Site), "NL", "FIN"))
  res_temp <- c()
  for(k in unique(dat.Col$Species)){
    dat.Col.temp <- dat.Col %>% filter(Species == k)
    dat.Ext.temp <- dat.Ext %>% filter(Species == k)
    
    for(l in unique(dat.Col.temp$Country)){
      res_temp2 <- cbind(
        Species = k, 
        Country = l,
        thresholds = paste(i, "-", j),
        dat.Col.temp %>% filter(Country == l) %>% ungroup() %>% 
          summarise(n.col = length(Colonisation[Colonisation == 1 & !is.na(Colonisation)]), 
                    n.NonCol = length(Colonisation[Colonisation == 0 & !is.na(Colonisation)]), 
                    n.indeterCol = length(Colonisation[is.na(Colonisation)]),
                    prob.col = n.col /(n.col + n.NonCol)),
        
        dat.Ext.temp %>% filter(Country == l) %>% ungroup() %>% 
          summarise(n.ext = length(Extinction[Extinction == 1 & !is.na(Extinction)]), 
                    n.NonExt = length(Extinction[Extinction == 0 & !is.na(Extinction)]), 
                    n.indeterExt = length(Extinction[is.na(Extinction)]),
                    prob.ext = n.ext /(n.ext + n.NonExt))
      )
      res_temp <- rbind(res_temp, res_temp2)
    }
  }
  return(res_temp)
}

table_ColExt_sp <- table_ColExt_sp %>% group_by(Country, thresholds) %>% mutate(rank.col = dense_rank(desc(prob.col)), 
                                                                                rank.ext = dense_rank(desc(prob.ext)))


# table
table_ColExt_sp2 <- c()
for(i in unique(table_ColExt_sp$Country)){
  res3 <- data.frame(Variables = unique(table_ColExt_sp[table_ColExt_sp$Country == i, "Species"]), Country = i)
  for(k in unique(table_ColExt_sp$thresholds)){
    res_temp <- table_ColExt_sp %>% filter(Country == i,  thresholds == k)
    colnames(res_temp)[4:13] <- paste(k, colnames(res_temp)[4:13])
    res3 <- cbind(res3, res_temp[,4:13])
  }
  table_ColExt_sp2 <- rbind(table_ColExt_sp2, res3)
}

table_ColExt_sp3 <- full_join(table_ColExt_sp2 %>% filter(Country == "NL"), 
                              table_ColExt_sp2 %>% filter(Country == "FIN"), 
                              by = "Species")


write.table(table_ColExt_sp3, "clipboard_124000")


registerDoFuture()
options(future.globals.maxSize = Inf)
plan(multiprocess, workers = 4)

res <- foreach(i = c(1:4)/10, j = rev(c(6:9)/10), .combine = rbind.data.frame) %:% foreach(k = unique(butterflies$Scale), .combine = rbind.data.frame) %dopar%{
  
  dat <- sp_site_occTrend(left_join(dat.occ.trend, butterflies %>% dplyr::select(2,8) %>% group_by(Site) %>% summarise_all(first)), i, j)
  dat <- left_join(dat, butterflies %>% dplyr::select(2,3,5,6,8,10:22) %>% group_by(Species, Site, Scale) %>% summarise_all(first))
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
res$Variables <- gsub("PC3", "Generation time", res$Variables)
res$Variables <- gsub("PC4", "Resource specialisation", res$Variables)
res$Variables <- gsub("PLAND", "Proportion of SNH", res$Variables)
res$Variables <- gsub("CLUMPY", "Aggregation of SNH", res$Variables)
res$Variables <- gsub(":", " × ", res$Variables)
res$Variables <- gsub("Proportion of SNH × Aggregation of SNH", "Proportion × Aggregation of SNH", res$Variables)
res$Variables <- gsub("X", "Longitude", res$Variables)
res$Variables <- gsub("Y", "Latitude", res$Variables)

res$Variables <- gsub("STI × Proportion × Aggregation of SNH", "STI × Proportion ×\nAggregation of SNH", res$Variables)
res$Variables <- gsub("Spatial use × Proportion × Aggregation of SNH", "Spatial use × Proportion ×\nAggregation of SNH", res$Variables)


res <- res %>% mutate(Significancy = ifelse(lwr < 0 & upr < 0 | lwr > 0 & upr > 0, "Significant", "Not significant"))

res$Scale <- factor(res$Scale, levels = rev(paste(c(1,3,5,10,20,30,50), "km")))

res$Variables <- factor(res$Variables, levels = c("(Intercept)",
                                                  "STI", "Spatial use", "Generation time", "Resource specialisation", 
                                                  "STI × Spatial use", "STI × Generation time", "STI × Resource specialisation",
                                                  "HabitatGeneralist", "HabitatOpen",
                                                  "Proportion of SNH", "Aggregation of SNH", "Proportion × Aggregation of SNH",
                                                  "Spatial use × Proportion of SNH", "Spatial use × Aggregation of SNH", 
                                                  "Spatial use × Proportion ×\nAggregation of SNH",
                                                  "STI × Proportion of SNH", "STI × Aggregation of SNH", 
                                                  "STI × Proportion ×\nAggregation of SNH", 
                                                  "Longitude", "Latitude", "Longitude × Latitude"))


res <- res %>% arrange(Variables, Scale, Process) %>% mutate(Group = rep(c(rep("Various", 1), 
                                                                           rep("Traits", 7),
                                                                           rep("Various", 2),
                                                                           rep("Fragmentation", 3),
                                                                           rep("Fragmentation × traits", 6),
                                                                           rep("Various", 3)), each = 56))

res$Group <- factor(res$Group, levels = c("Various", "Traits", "Fragmentation", "Fragmentation × traits"))

### plot ###

ggplot(res %>% dplyr::filter(grepl("STI|Spatial use|Generation time|Resource specialisation|Proportion|Aggregation", Variables), 
                             Range == "0.2 - 0.8") %>% 
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

# table
res2 <- c()
for(i in unique(res$Process)){
  res4 <- c()
  for(j in unique(res$Range)){
    res3 <- data.frame(Variables = unique(res$Variables), Process = i, Range = j)
    for(k in unique(res$Scale)){
      res_temp <- res %>% filter(Process == i, Range == j, Scale == k)
      colnames(res_temp)[5:7] <- paste(k, colnames(res_temp)[5:7])
      res3 <- cbind(res3, res_temp[,5:7])
    }
    res4 <- rbind(res4, res3)
  }
  res2 <- rbind(res2, res4)
}

write.table(res2 %>% arrange(Process, Range, Variables), "clipboard-16384", row.names = F)

###############################
### drivers of colonisation ###
###     and extinction at   ###
###      20km scale and     ###
###    .2 - .8 threshold    ###
###############################


### extract data at scale = 20km and range = 0.2 - 0.8 ###
dat <- sp_site_occTrend(left_join(dat.occ.trend, butterflies %>% dplyr::select(2,8) %>% group_by(Site) %>% summarise_all(first)), .2, .8)
table(dat$Trend)

dat <- left_join(dat, butterflies %>% dplyr::select(2,3,5,6,8,10:22) %>% 
                   group_by(Species, Site, Scale) %>% summarise_all(first))
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
library(DHARMa)

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

testDispersion(simulateResiduals(m.col))

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

testDispersion(simulateResiduals(m.ext))

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
       xtrans = function(x){x * scaleList$scale["PC1"] + scaleList$center["PC1"]}, ylim = c(0,1),
       xlab = "Spatial use", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

visreg(m.col, xvar = "PC4", by  ="STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = T,
       xtrans = function(x){x * scaleList$scale["PC4"] + scaleList$center["PC4"]}, ylim = c(0,1),
       xlab = "Resource specialisation", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

visreg(m.col, xvar = "PC3", by  ="STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = T,
       xtrans = function(x){x * scaleList$scale["PC3"] + scaleList$center["PC3"]}, ylim = c(0,1),
       xlab = "Generation time", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

# plot(0,0,axes = F, ann = F, col = "white")
# legend(x = -1, y = 0, title = "STI", legend = c("Low", "High"), lty = 1, col = c("#1C86EE", "#EEC900"), cex = 1.2, lwd = 1.2, bty = "n")

#ext

visreg(m.ext, xvar = "PC1", by = "STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = T,
       xtrans = function(x){x * scaleList$scale["PC1"] + scaleList$center["PC1"]},ylim = c(0,1),
       xlab = "Spatial use", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

visreg(m.ext, xvar = "PC3", by  ="STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = T,
       xtrans = function(x){x * scaleList$scale["PC3"] + scaleList$center["PC3"]}, ylim = c(0,1),
       xlab = "Resource specialisation", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

visreg(m.ext, xvar = "PC4", by  ="STI_rel", scale = "response", rug = F, breaks = c(-1,1), overlay = T, band = T,
       xtrans = function(x){x * scaleList$scale["PC4"] + scaleList$center["PC4"]}, ylim = c(0,1),
       xlab = "Generation time", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

dev.off()

# landscape x traits

# 1
l1 <- predict_raster2(m.col, xvar = "PLAND", yvar = "CLUMPY", scaleList = scaleList, cond = list(STI_rel = -1))
l2 <- predict_raster2(m.ext, xvar = "PLAND", yvar = "CLUMPY", scaleList = scaleList, cond = list(STI_rel = -1))

l3 <- predict_raster2(m.col, xvar = "PLAND", yvar = "CLUMPY", scaleList = scaleList, cond = list(STI_rel = 1))
l4 <- predict_raster2(m.ext, xvar = "PLAND", yvar = "CLUMPY", scaleList = scaleList, cond = list(STI_rel = 1))

landTraitEffect <- bind_rows(
  "Low STI" = bind_rows(Colonisation = l1, Extinction = l2, .id = "Process"),
  "High STI" = bind_rows(Colonisation = l3, Extinction = l4, .id = "Process"),
  .id = "STI")

p1 <- ggplot(landTraitEffect %>% filter(Process == "Colonisation"), aes(x = PLAND, y = CLUMPY, fill = pred, z = pred)) +
  geom_tile() + facet_grid(STI~., scales = "free", drop = F) +
  scale_fill_gradient(name = "Colonisation\nprobability\n",
                      low=c("gold2"), high = "dodgerblue") +
  scale_x_continuous("Proportion of SNH") + scale_y_continuous("Aggregation of SNH") +
  stat_contour(color="white", size=0.1) +
  theme_classic() +
  theme(panel.border=element_blank(),
        strip.text=element_text(size=12, colour="black"),
        strip.background=element_rect(colour="white", 
                                      fill="white"))

p2 <- ggplot(landTraitEffect %>% filter(Process == "Extinction"), aes(x = PLAND, y = CLUMPY, fill = pred, z = pred)) +
  geom_tile() + facet_grid(STI~., scales = "free", drop = F) +
  scale_fill_gradient(name = "Extinction\nprobability\n",
                      high=c("gold2"), low = "dodgerblue") +
  scale_x_continuous("Proportion of SNH") + scale_y_continuous("Aggregation of SNH") +
  stat_contour(color="white", size=0.1) +
  theme_classic() +
  theme(panel.border=element_blank(),
        strip.text=element_text(size=12, colour="black"),
        strip.background=element_rect(colour="white", 
                                      fill="white"))

p3 <- cowplot::plot_grid(p1, p2, labels = c("Colonisation", "Extinction"))

cowplot::ggsave("../plot_traitsxLand3.pdf", p3, width = 8, height = 4.5, dpi = 300)

# 2
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

p1 <- ggplot(landTraitEffect %>% filter(Process == "Colonisation"), aes(x= STI_rel, y = PC1, fill = estimate)) +
  geom_tile() + facet_grid(Land_var~., scales = "free", drop = F) +
  scale_fill_gradient2(name = "Effect on \ncolonisation\nprobability\n",
                       low=c("gold2"), high = "dodgerblue", mid  ="light grey") +
  scale_x_continuous("Species temperature index") + scale_y_continuous("Spatial use") +
  theme_classic() +
  theme(panel.border=element_blank(),
        strip.text=element_text(size=12, colour="black"),
        strip.background=element_rect(colour="white", 
                                      fill="white"))

p2 <- ggplot(landTraitEffect %>% filter(Process == "Extinction"), aes(x= STI_rel, y = PC1, fill = estimate)) +
  geom_tile() + facet_grid(Land_var~., scales = "free") +
  scale_fill_gradient2(name = "Effect on \nextinction\nprobability\n",
                       low=c("dodgerblue"), high = "gold2", mid  ="light grey") +
  scale_x_continuous("Species temperature index") + scale_y_continuous("Spatial use") +
  theme_classic() +
  theme(panel.border=element_blank(),
        strip.text=element_text(size=12, colour="black"),
        strip.background=element_rect(colour="white", 
                                      fill="white"))

p3 <- cowplot::plot_grid(p1, p2, labels = c("Colonisation", "Extinction"))

cowplot::ggsave("../plot_traitsxLand3.pdf", p3, width = 8, height = 4.5, dpi = 300)


# alt [1]
par(mfrow=c(2,2))
visreg(m.col, xvar = "PLAND", by = "STI_rel",
       breaks = seq(min(range(dat.Col$STI_rel, na.rm = T)),
                    max(dat.Col$STI_rel, na.rm = T), length.out = 50), overlay = T, 
       scale = "response", rug = F, ylim = c(0,1), band = F, legend = F,
       line.par=list(col=colorRampPalette(c("#1C86EE", "#EEC900"))(50)))
visreg(m.col, xvar = "CLUMPY", by = "STI_rel",
       breaks = seq(min(range(dat.Col$STI_rel, na.rm = T)),
                    max(dat.Col$STI_rel, na.rm = T), length.out = 50), overlay = T, 
       scale = "response", rug = F, ylim = c(0,1), band = F, legend = F,
       line.par=list(col=colorRampPalette(c("#1C86EE", "#EEC900"))(50)))

visreg(m.ext, xvar = "PLAND", by = "STI_rel",
       breaks = seq(min(range(dat.Col$STI_rel, na.rm = T)),
                    max(dat.Col$STI_rel, na.rm = T), length.out = 50), overlay = T, 
       scale = "response", rug = F, ylim = c(0,1), band = F, legend = F,
       line.par=list(col=colorRampPalette(c("#1C86EE", "#EEC900"))(50)))
visreg(m.ext, xvar = "CLUMPY", by = "STI_rel",
       breaks = seq(min(range(dat.Col$STI_rel, na.rm = T)),
                    max(dat.Col$STI_rel, na.rm = T), length.out = 50), overlay = T, 
       scale = "response", rug = F, ylim = c(0,1), band = F, legend = F,
       line.par=list(col=colorRampPalette(c("#1C86EE", "#EEC900"))(50)))
par(mfrow=c(1,1))

# alt [2]
par(mfrow=c(2,2))
visreg(m.col, xvar = "PLAND", by = "STI_rel",
       breaks = c(-1,1), overlay = T, 
       scale = "response", rug = F, ylim = c(0,1), band = T,
       xtrans = function(x){x * scaleList$scale["PLAND"] + scaleList$center["PLAND"]},ylim = c(0,1),
       xlab = "Proportion of SNH", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

visreg(m.col, xvar = "CLUMPY", by = "STI_rel",
       breaks = c(-1,1), overlay = T, 
       scale = "response", rug = F, ylim = c(0,1), band = T,
       xtrans = function(x){x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"]},ylim = c(0,1),
       xlab = "Aggregation of SNH", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

visreg(m.ext, xvar = "PLAND", by = "STI_rel",
       breaks = c(-1,1), overlay = T, 
       scale = "response", rug = F, ylim = c(0,1), band = T,
       xtrans = function(x){x * scaleList$scale["PLAND"] + scaleList$center["PLAND"]},ylim = c(0,1),
       xlab = "Proportion of SNH", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

visreg(m.ext, xvar = "CLUMPY", by = "STI_rel",
       breaks = c(-1,1), overlay = T, 
       scale = "response", rug = F, ylim = c(0,1), band = T,
       xtrans = function(x){x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"]},ylim = c(0,1),
       xlab = "Aggregation of SNH", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))
par(mfrow=c(1,1))

# alt [3] spatial use
par(mfrow=c(2,2))
visreg(m.col, xvar = "PLAND", by = "PC1",
       breaks = c(-1,1), overlay = T, 
       scale = "response", rug = F, ylim = c(0,1), band = T,
       xtrans = function(x){x * scaleList$scale["PLAND"] + scaleList$center["PLAND"]},ylim = c(0,1),
       xlab = "Proportion of SNH", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

visreg(m.col, xvar = "CLUMPY", by = "PC1",
       breaks = c(-1,1), overlay = T, 
       scale = "response", rug = F, ylim = c(0,1), band = T,
       xtrans = function(x){x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"]},ylim = c(0,1),
       xlab = "Aggregation of SNH", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

visreg(m.ext, xvar = "PLAND", by = "PC1",
       breaks = c(-1,1), overlay = T, 
       scale = "response", rug = F, ylim = c(0,1), band = T,
       xtrans = function(x){x * scaleList$scale["PLAND"] + scaleList$center["PLAND"]},ylim = c(0,1),
       xlab = "Proportion of SNH", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))

visreg(m.ext, xvar = "CLUMPY", by = "PC1",
       breaks = c(-1,1), overlay = T, 
       scale = "response", rug = F, ylim = c(0,1), band = T,
       xtrans = function(x){x * scaleList$scale["CLUMPY"] + scaleList$center["CLUMPY"]},ylim = c(0,1),
       xlab = "Aggregation of SNH", legend = F, fill.par=list(col=c("#1C86EE60", "#EEC90060")),
       line.par=list(col=c("#1C86EE", "#EEC900")))
par(mfrow=c(1,1))




###########################
## remove unused objects ##
###########################


gdata:::keep(butterflies, dat.occ.trend,
             sure = T)

