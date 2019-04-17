library(data.table)
library(tidyverse)
library(raster)
library(sp)
library(landscapemetrics)
library(rgeos)
library(randomForest)
library(cooccur)

#####################
## background info ##
#####################

expanding_species <- countsFIN %>% group_by(Year) %>% mutate(nSites_tot = length(unique(Site))) %>% group_by(Species, Year) %>%
  summarise(nSites = length(unique(Site)), nSites_tot = unique(nSites_tot)) %>% group_by(Species) %>%
  do(trend = glm(cbind(nSites,nSites_tot) ~ Year, family = binomial, data = .)) %>% 
  broom::tidy(trend) %>% filter(term == "Year") %>% filter(p.value < .05, estimate > 0) %>% arrange(desc(estimate)) %>%
  pull(Species)


#################
### load data ###
#################

##### butterfly data ######
countsFIN1 <- fread("../Data/Butterflies - Finland/FINLAND_Records_1999-2015.txt", sep = ";", h=T)
countsFIN2 <- fread("../Data/Butterflies - Finland/FINLAND_Records_2016.txt", sep = ";", h=T)[,-6]
countsFIN <- rbind(countsFIN1, countsFIN2)
colnames(countsFIN) <- gsub ("Species_Faunaeur", "Species", colnames(countsFIN))
colnames(countsFIN) <- gsub ("Individuals", "n", colnames(countsFIN))

countsFIN <- countsFIN %>% 
  group_by(Site, Species, Year) %>% 
  summarise(n = max(n)) %>% as.tbl # keep max no. individuals

#############################
##### extract site data #####
#############################
FIN_wgs84 <- getData("GADM", country = "FIN", level = 0)
FIN_3035 <- spTransform(FIN_wgs84, crs("+init=epsg:3035"))

## load land use data ##
lc_class <- as.tbl(rio::import(file = "../Landcover/Corine_land-cover_2012_raster/clc_legend.xls", which = 1L))

LU <- raster("../Landcover/Corine_land-cover_2012_raster/g100_clc12_V18_5.tif")
LU <- crop(LU, extent(FIN_3035)*1.1)

## load climate data ##
pre <- stack(lapply(list.files("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/CRU_data\\pre", full.names = T), function(x)(stack(x))))
tmp <- stack(lapply(list.files("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/CRU_data\\tmp", full.names = T), function(x)(stack(x))))
tmn <- stack(lapply(list.files("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/CRU_data\\tmn", full.names = T), function(x)(stack(x))))
tmx <- stack(lapply(list.files("//storage.slu.se/Home$/yofo0001/My Documents/Recherche/CRU_data\\tmx", full.names = T), function(x)(stack(x))))

## site data ##
sites <- read_csv(file = "../Data/Butterflies - Finland/Sites_FIN_ETRS89_landcover.csv")
sites_3035 <- SpatialPoints(sites[,c("X", "Y")], proj4string=CRS("+init=epsg:3035"))
sites_wgs84 <- spTransform(sites_3035, crs(clim))

## background 
# back_3035 <- spsample(FIN_3035, 1000, type = "random")
# back_wgs84 <- spTransform(back_3035, crs(LU))

## extract temperature data
sites_temp <- c()
for(z in unique(countsFIN$Year)){
  
  pre_site_temp <- extract(mean(subset(pre, grep(z, names(pre)))), sites_wgs84)
  tmp_site_temp <- extract(mean(subset(tmp, grep(z, names(tmp)))), sites_wgs84)
  tmn_site_temp <- extract(mean(subset(tmn, grep(z, names(tmn)))), sites_wgs84)
  tmx_site_temp <- extract(mean(subset(tmx, grep(z, names(tmx)))), sites_wgs84)
  
  sites_temp <- rbind.data.frame(
    cbind.data.frame(sites[,1], Year = z,
                     pre = pre_site_temp, tmp = tmp_site_temp, tmn = tmn_site_temp, tmx = tmx_site_temp)
  )
}

# back_temp <- extract(clim, back_wgs84)

## extract land use data
# sites_LU <- extract(LU, sites_3035)
# 
# 
# sites_LU <- c()
# for(j in 1:nrow(sites)){
#   
#   for(k in c(1000,5000,10000)){
#     
#     buf <- gBuffer(sites_3035[j,], width = k)
#     land <- crop(LU, buf)
#     land <- mask(land, buf)
#     
#     sites_LU <- bind_rows(sites_LU,
#                           cbind.data.frame(sites[j,], Scale = k,
#                                            land %>% lsm_p_area() %>% mutate(tot = sum(value)) %>%
#                                              group_by(class) %>% summarise(value = sum(value)/unique(tot))))
#     
#   }
# }
# sites_LU <- sites_LU %>% left_join(lc_class, by = c("class" = "GRID_CODE")) %>% as.tbl %>%
#   group_by(Site, Scale, LABEL1) %>% summarise(landcover_1km = sum(value)) %>% spread("Scale", "landcover_1km", fill = 0)
# 
# colnames(sites_LU)[3:5] <- paste0("Landcover_", colnames(sites_LU)[3:5])
# 
# sites_LU_1000 <- sites_LU  %>% dplyr::select(c(1:3)) %>% spread("LABEL1", "Landcover_1000", fill = 0)
# colnames(sites_LU_1000)[-1] <- paste0(colnames(sites_LU_1000)[-1], "_1000")
# colnames(sites_LU_1000) <- gsub(" ","_", colnames(sites_LU_1000))
# 
# back_LU <- c()
# for(j in 1:nrow(back_3035@coords)){
#   
#   for(k in c(1000,5000,10000)){
#     
#     buf <- gBuffer(back_3035[j,], width = k)
#     land <- crop(LU, buf)
#     land <- mask(land, buf)
#     
#     
#     back_LU <- bind_rows(back_LU,
#                          cbind.data.frame(t(back_3035@coords[j,]), Scale = k,
#                                           land %>% lsm_p_area() %>% mutate(tot = sum(value)) %>%
#                                             group_by(class) %>% summarise(value = sum(value)/unique(tot))))
#     
#   }
# }
# back_LU <- back_LU %>% left_join(lc_class, by = c("class" = "GRID_CODE")) %>% as.tbl %>%
#   group_by(coords = paste(x, y), Scale, LABEL1) %>% summarise(landcover_1km = sum(value)) %>% 
#   spread("Scale", "landcover_1km", fill = 0) 
# 
# colnames(back_LU)[3:5] <- paste0("Landcover_", colnames(back_LU)[3:5])
# 
# back_LU_1000 <- back_LU  %>% dplyr::select(c(1:3)) %>% spread("LABEL1", "Landcover_1000", fill = 0)
# colnames(back_LU_1000)[-1] <- paste0(colnames(back_LU_1000)[-1], "_1000")
# colnames(back_LU_1000) <- gsub(" ","_", colnames(back_LU_1000))
# 
# 
# sites_LU_5000 <- sites_LU  %>% dplyr::select(c(1,2,4)) %>% spread("LABEL1", "Landcover_5000", fill = 0)
# colnames(sites_LU_5000)[-1] <- paste0(colnames(sites_LU_5000)[-1], "_5000")
# 
# sites_LU_10000 <- sites_LU %>% dplyr::select(c(1,2,5)) %>% spread("LABEL1", "Landcover_10000", fill = 0)
# colnames(sites_LU_10000)[-1] <- paste0(colnames(sites_LU_10000)[-1], "_10000")
# 
# sites_LU <- bind_cols(sites_LU_1000, sites_LU_5000[,-1]) %>% bind_cols(sites_LU_10000[,-1])
# colnames(sites_LU) <- gsub(" ","_", colnames(sites_LU))

vars <- c("bio1", "bio4", "bio5", "bio6", "bio12", "bio13", "bio14", "bio15")

sites_env <- cbind.data.frame(sites_LU_1000, sites_temp[,vars])
sites_clim_pca <- cbind.data.frame(sites[,1], ade4::dudi.pca(sites_env[,-c(1:7)], scannf = F)$li) %>% left_join(sites)

# back_env <- cbind.data.frame(back_LU_1000[,-1], back_temp[,vars])

###########################################
### compute species- and site-specific  ###
###         habitat suitability         ###
###########################################


# sp_site_suitability <- c()
# for(i in unique(countsFIN$Species)){
#   
#   data_temp_pres <- countsFIN %>% ungroup %>% filter(Species == i) %>% left_join(sites) %>% dplyr::select(1,5,6) %>%
#     group_by(Site) %>% summarise_all(unique)
#   # data_temp_abs <- countsFIN %>% filter(!Site %in% data_temp_pres$Site) %>% group_by(Site) %>% summarise(occ = 0)
#   # data_temp <- bind_rows(data_temp_pres,data_temp_abs)
#   
#   # m <- maxent(x = sites_env %>% filter(Site %in% data_temp_pres$Site) %>% dplyr::select(-1) %>% 
#   #               bind_rows(back_env),
#   #             p = c(rep(1, nrow(data_temp_pres)), rep(0, nrow(back_env))))
#   
#   d <- cbind(occ = c(rep(1, nrow(data_temp_pres)), rep(0, nrow(back_env))), 
#              sites_env %>% filter(Site %in% data_temp_pres$Site) %>% dplyr::select(-1) %>% 
#                bind_rows(back_env))
#   
#   m <- randomForest(as.factor(occ) ~ ., data = d,
#                     sampsize = rep(min(table(d$occ)), 2))
#   
#   sp_site_suitability <- rbind.data.frame(
#     sp_site_suitability,
#     cbind.data.frame(
#       Species = i, 
#       Site = sites_env[,1],
#       suitable = c(rep(1, nrow(data_temp_pres)),
#                    as.numeric(as.character(predict(m, 
#                                                    newdata = sites_env %>% 
#                                                      filter(!Site %in% data_temp_pres$Site)))))
#     )
#   )
#   
# }


#######################################
### compute non-random associations ###
#######################################

assoc.Sp2 <- c()
for(m in c(min(countsFIN$Year) : max(countsFIN$Year))){
  countsFIN_comm <- countsFIN %>% ungroup() %>% filter(Year == m) %>% mutate(n = 1) %>% 
    spread(Species, n, fill = 0) %>% as.data.frame %>% 
    mutate(id = paste0(Site,"_",Year)) %>% dplyr::select(-1,-2) %>% column_to_rownames("id")
  countsFIN_comm <- countsFIN_comm[,colSums(countsFIN_comm) > (0.01 * nrow(countsFIN_comm)) & 
                                     colSums(countsFIN_comm) < (0.99 * nrow(countsFIN_comm))]
  
  assoc.Sp_temp <- cooccur(t(countsFIN_comm), type = "spp_site", thresh = F, spp_names = T,
                           prob = "comb")
  assoc.Sp_temp <- assoc.Sp_temp$results
  assoc.Sp_temp <- assoc.Sp_temp %>% mutate(association = ifelse(p_gt < 0.05, "Aggregated", 
                                                                 ifelse(p_lt < 0.05, "Segregated", "Random")),
                                            Year = m)
  assoc.Sp2 <- rbind.data.frame(assoc.Sp2, assoc.Sp_temp)
}
assoc.Sp2 <- as.tbl(assoc.Sp2)

# assoc.Sp3 <- assoc.Sp2 %>% mutate(pair = paste(sp1_name, "-", sp2_name)) %>% group_by(pair, association) %>% 
#   mutate(n.association = n()) %>% group_by(pair) %>% mutate(n.yr = n()) %>%
#   summarise(association = association[n.association = max(n.association)], 
#             prop.association = max(n.association)/unique(n.yr),
#             n.yr = unique(n.yr)) %>%
#   mutate(association = ifelse(prop.association < 0.7, "Random", association)) %>% filter(n.yr > 4) %>%
#   separate(pair, c("sp1", "sp2"), sep = " - ")

table(assoc.Sp3$association)

## with species-site mask ##

sp_site_suitability2 <- countsFIN %>% dplyr::select(-4) %>% ungroup %>% 
  complete(Species, Year, Site) %>%
  left_join(sp_site_suitability) %>% mutate(suitable = ifelse(suitable == 1, 1, 0)) %>% 
  filter(Species %in% colnames(countsFIN_comm)) %>%
  group_by(Species, Site) %>% spread(Species, suitable) %>% as.data.frame %>% 
  mutate(id = paste0(Site,"_",Year)) %>% dplyr::select(-1,-2)%>% 
  filter(id %in% rownames(countsFIN_comm))  %>% column_to_rownames("id")

assoc.Sp <- cooccur(t(countsFIN_comm), type = "spp_site", thresh = F, spp_names = T,
                    site_mask = t(sp_site_suitability2),
                    prob = "comb")
plot(assoc.Sp)
pair.profile(assoc.Sp)
obs.v.exp(assoc.Sp)


## filter interactions ##

res_assoc <- c()
for(l in 1:nrow(assoc.Sp2)){
  
  diff_env <- NULL
  diff_dist <- NULL
  
  sp1 <- c(assoc.Sp2 %>% slice(l) %>% pull(sp1_name) %>% as.character)
  sp2 <- c(assoc.Sp2 %>% slice(l) %>% pull(sp2_name) %>% as.character)
  
  
  if(assoc.Sp2 %>% slice(l) %>% pull(association) == "Segregated"){
    
    assoc_cond_10 <- countsFIN %>% ungroup() %>% filter(Species %in% sp1, 
                                                        Year == assoc.Sp2 %>% slice(l) %>% pull(Year)) %>% 
      filter(!Species %in% sp2) %>% 
      dplyr::select(Site, Year) %>% distinct %>%
      left_join(sites_env_pca)
    
    assoc_cond_01 <- countsFIN %>% ungroup() %>% filter(Species %in% sp2, 
                                                        Year == assoc.Sp2 %>% slice(l) %>% pull(Year)) %>% 
      filter(!Species %in% sp1) %>% 
      dplyr::select(Site, Year) %>% distinct %>%
      left_join(sites_env_pca)
    
    assoc_cond_merged <- bind_rows("10" = assoc_cond_10, "00" = assoc_cond_01, .id = "assoc")
    
  } else if(assoc.Sp2 %>% slice(l) %>% pull(association) == "Aggregated"){
    
    assoc_cond_11 <- countsFIN %>% ungroup() %>% filter(Species %in% c(sp1, sp2), 
                                                        Year == assoc.Sp2 %>% slice(l) %>% pull(Year)) %>% 
      group_by(Site, Year) %>% summarise(n = length(unique(Species))) %>% filter(n > 1) %>%
      dplyr::select(Site, Year) %>% distinct %>%
      dplyr::select(1) %>% 
      left_join(sites_env_pca)
    
    assoc_cond_00 <- countsFIN %>% ungroup() %>% filter(!Species %in% c(sp1, sp2), 
                                                        !Site %in% assoc_cond_11$Site, 
                                                        Year == assoc.Sp2 %>% slice(l) %>% pull(Year)) %>% 
      dplyr::select(Site, Year) %>% distinct %>%
      left_join(sites_env_pca)
    
    assoc_cond_merged <- bind_rows("00" = assoc_cond_00, "11" = assoc_cond_11, .id = "assoc")
    
  }
  
  association <- assoc.Sp2 %>% slice(l) %>% pull(association)
  
  if(association != "Random"){
    diff_env <- summary(manova(cbind(Axis1, Axis2) ~ assoc, data = assoc_cond_merged))$stat[1,6]
    diff_dist <- summary(manova(cbind(X, Y) ~ assoc, data = assoc_cond_merged))$stat[1,6]
  }
  
  
  if(association == "Random"){
    explanation <- "Random"
  } else if(association == "Aggregated"){
    if(diff_dist < 0.05){
      if(diff_env < 0.05){
        explanation <- "Dispersal_limitation_or_environmental_filtering"
      } else {
        explanation <- "Dispersal_limitation"
      }
    } else {
      if(diff_env < 0.05){
        explanation <- "Environmental_filtering"
      } else {
        explanation <- "Positive_species_interaction"
      }
    }
  } else if(association == "Segregated"){
    if(diff_dist < 0.05){
      if(diff_env < 0.05){
        explanation <- "Dispersal_limitation_or_environmental_filtering"
      } else {
        explanation <- "Dispersal_limitation"
      }
    } else {
      if(diff_env < 0.05){
        explanation <- "Environmental_filtering"
      } else {
        explanation <- "Negative_species_interaction"
      }
    }
  }
  
  res_assoc <- rbind.data.frame(res_assoc, 
                                cbind.data.frame(sp1 = sp1, sp2 = sp2,
                                                 Year = assoc.Sp2[l,"Year"],
                                                 association, explanation)
  )
}

res_assoc <- as.tbl(data.frame(res_assoc))
table(res_assoc$explanation)

res_assoc[res_assoc$explanation == "Negative_species_interaction",]
res_assoc[res_assoc$explanation == "Positive_species_interaction",]

############################
## change in cooccurrence ##
############################
library(emmeans)
library(lmerTest)

countsFIN %>% group_by(Site, Year) %>% summarise(l = list(unique(Species)))

sum.assoc1 <- res_assoc %>% left_join(countsFIN, by = c("sp1" = "Species", "Year" = "Year")) %>% 
  group_by(Site, Year, explanation) %>% 
  summarise(n.explanation = n())
sum.assoc2 <- res_assoc %>% left_join(countsFIN, by = c("sp2" = "Species", "Year" = "Year")) %>% 
  group_by(Site, Year, explanation) %>% 
  summarise(n.explanation = n())
sum.assoc <- bind_rows(sum.assoc1, sum.assoc2) %>% distinct()

sum.assoc.yr <- sum.assoc %>% group_by(Year) %>% mutate(n.sites = length(unique(Site))) %>%
  group_by(Year, explanation) %>% summarise(n = sum(n.explanation)/unique(n.sites))

sum.assoc.yr %>% filter(explanation %in% c("Positive_species_interaction", "Negative_species_interaction")) %>%
  ggplot(aes(x = Year, y = n, color = explanation)) + geom_point() + geom_line() + 
  geom_smooth() +
  facet_wrap(~explanation, scales = "free")

sum.assoc.site <- sum.assoc %>% left_join(sites) %>% left_join(sites_env)

## temporal trend ##
m.pos <- glmer(n.explanation ~ Year + (1|Site), family = gaussian(link = log),
               data = sum.assoc %>% filter(explanation == "Positive_species_interaction") %>% 
                 ungroup %>% mutate(Year = Year - min(Year)))
summary(m.pos)

m.neg <- glmer(n.explanation ~ Year + (1|Site), family = gaussian(link = log),
               data = sum.assoc %>% filter(explanation == "Negative_species_interaction") %>% 
                 ungroup %>% mutate(Year = Year - min(Year)))
summary(m.neg)

## spatial * temporal trend ##
m.pos.geo <- glmer(n.explanation ~ Year * (X*Y) + (1|Site), family = gaussian(link = log),
                   data = sum.assoc.site %>% filter(explanation == "Positive_species_interaction") %>% ungroup %>%
                     mutate(X = scale(X), Y = scale(Y), Year = Year - min(Year)))
summary(m.pos.geo)

m.neg.geo <- glmer(n.explanation ~ Year * (X*Y) + (1|Site), family = gaussian(link = log),
                   data = sum.assoc.site %>% filter(explanation == "Negative_species_interaction") %>% ungroup %>%
                     mutate(X = scale(X), Y = scale(Y), Year = Year - min(Year)))
summary(m.neg.geo)

## effect of climate ##
dat_temp <- sum.assoc.site %>% filter(explanation == "Positive_species_interaction") %>% ungroup %>%
  mutate(X = scale(X), Y = scale(Y))
m.pos.env <- glmer(n.explanation ~ X*Y + 
                     bio1 + bio4 + bio12 + bio15 + 
                     (1|Site), family = gaussian(link = log),
                   data = dat_temp, na.action = na.fail)
summary(m.pos.env)
d.pos.env <- dredge(m.pos.env)
avg.pos.env <- model.avg(d.pos.env, subset = delta < 2, fit = T)


dat_temp2 <- sum.assoc.site %>% filter(explanation == "Negative_species_interaction") %>% ungroup %>%
  mutate(X = scale(X), Y = scale(Y))
m.neg.env <- glmer(n.explanation ~ X*Y + 
                     bio1 + bio4 + bio12 + bio15 + 
                     (1|Site), family = gaussian(link = log),
                   data = dat_temp2, na.action = na.fail)
summary(m.neg.env)
d.neg.env <- dredge(m.neg.env)
avg.neg.env <- model.avg(d.neg.env, subset = delta < 2, fit = T)



## plots ##
# spatial trend
FIN_ras <- projectRaster(clim[[1]], crs = crs(sites_3035))
FIN_ras <- mask(FIN_ras, FIN_3035)

FIN_ras_x <- FIN_ras
FIN_ras_x <- setValues(FIN_ras_x, xyFromCell(FIN_ras, 1:ncell(FIN_ras))[,1])
names(FIN_ras_x) <- "X"

FIN_ras_y <- FIN_ras
FIN_ras_y <- setValues(FIN_ras_y, xyFromCell(FIN_ras, 1:ncell(FIN_ras))[,2])
names(FIN_ras_y) <- "Y"

FIN_ras <- stack(FIN_ras_x, FIN_ras_y)
FIN_ras <- mask(FIN_ras, FIN_3035)
FIN_ras <- scale(FIN_ras)


fun.trend=function(x) {if(is.na(x[1])){ NA } else {m = glm(x ~ c(1:18), family = gaussian(link=log));summary(m)$coefficients[2]}}

pred_time_series_pos <- lapply(0:17, function(x)predict(FIN_ras, m.pos.geo, const=(data.frame(Site=2, Year = x))))
slope_time_series_pos <- calc(stack(pred_time_series_pos), fun = fun.trend)
slope_time_series_pos_cropped <- crop(slope_time_series_pos, extent(sites_3035)*1.05)

pred_time_series_neg <- lapply(0:17, function(x)predict(FIN_ras, m.neg.geo, const=(data.frame(Site=2, Year = x))))
slope_time_series_neg <- calc(stack(pred_time_series_neg), fun = fun.trend)
slope_time_series_neg_cropped <- crop(slope_time_series_neg, extent(sites_3035)*1.05)


# spplot(
#   stack((exp(slope_time_series_pos_cropped) - 1) * 100,
#         (exp(slope_time_series_neg_cropped) - 1) * 100),
#   sp.layout = list("sp.points", sites_3035, pch=16, col='black'),
#   names.attr=c('Positive associations', 'Negative associations'),
#   col.regions = colorRampPalette(c("red", "blue")),
#   at = seq(-2, 0.0225, length.out = 100),
#   contour = T
# )


# Mean
par(mfrow=c(1,2))

plot(crop(predict(FIN_ras, m.pos.geo, const=(data.frame(Site=2, Year = 1))), extent(sites_3035)*1.05), 
     legend = T, axes = F, box = F,
     main = "Mean no. positive associations / site", cex.main = .8)
points(sites_3035@coords)

plot(crop(predict(FIN_ras, m.neg.geo, const=(data.frame(Site=2, Year = 1))), extent(sites_3035)*1.05), 
     legend = T, axes = F, box = F,
     main = "Mean no. negative associations / site", cex.main = .8)
points(sites_3035@coords)

par(mfrow=c(1,1))

# Trend
par(mfrow=c(1,2))

plot((exp(slope_time_series_pos_cropped) - 1) * 100, legend = T, axes = F, box = F,
     main = "% change in no. positive \nassociations / site / year", cex.main = .8)
points(sites_3035@coords)

plot((exp(slope_time_series_neg_cropped) - 1) * 100, legend = T, axes = F, box = F,
     main = "% change in no. positive \nassociations / site / year", cex.main = .8)
points(sites_3035@coords)

par(mfrow=c(1,1))


# temporal trends
m.pos.fac <- glmer(n.explanation ~ as.factor(Year) + (1|Site), family = gaussian(link = log),
                   data = sum.assoc %>% filter(explanation == "Positive_species_interaction"))
m.pos.cont <- glmer(n.explanation ~ Year + (1|Site), family = gaussian(link = log), 
                    data = sum.assoc %>% filter(explanation == "Positive_species_interaction"))
pred.pos.fac <- as.data.frame(emmeans(m.pos.fac, ~ Year))

m.neg.fac <- glmer(n.explanation ~ as.factor(Year) + (1|Site), family = gaussian(link = log),
                   data = sum.assoc %>% filter(explanation == "Negative_species_interaction"))
m.neg.cont <- glmer(n.explanation ~ Year + (1|Site), family = gaussian(link = log),
                    data = sum.assoc %>% filter(explanation == "Negative_species_interaction"))
pred.neg.fac <- as.data.frame(emmeans(m.neg.fac, ~ Year))


par(mfrow=c(1,3))

visreg(m.pos.cont, xvar = "Year", ylab = "No. positive associations / site", rug = F, scale = "response",
       ylim = range(exp(pred.pos.fac$emmean)))
points(exp(emmean) ~ Year, data = pred.pos.fac, type = "l")
points(exp(emmean) ~ Year, data = pred.pos.fac, pch = 16)

visreg(m.neg.cont, xvar = "Year", ylab = "No. negative associations / site", rug = F, scale = "response",
       ylim = range(exp(pred.neg.fac$emmean)))
points(exp(emmean) ~ Year, data = pred.neg.fac, type = "l")
points(exp(emmean) ~ Year, data = pred.neg.fac, pch = 16)

plot(pred.neg.fac$emmean, pred.pos.fac$emmean, ylab = "No. positive associations / site", 
     xlab = "No. negative associations / site")
abline(lm(pred.pos.fac$emmean~pred.neg.fac$emmean))

par(mfrow=c(1,1))


